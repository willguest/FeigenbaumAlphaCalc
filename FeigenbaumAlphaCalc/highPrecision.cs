using System;
using System.Collections.Generic;
using System.Linq;

using Extreme.Mathematics;
using Extreme.Mathematics.Calculus;
using Extreme.Mathematics.Generic;
using Extreme.Mathematics.LinearAlgebra;
using MathNet.Numerics.Interpolation;

namespace FeigenbaumAlphaCalc
{
    class highPrecision
    {
        int N = 700;
        int M = 1000000;
        BigFloat myalpha = 0;
        int count = 10;

        bool compute_alpha()
        {
            double[] betan = new double[N];
            betan[0] = -1.52763883147;
            betan[1] = 1.048327004372e-1;
            betan[2] = 2.669121419134012e-2;

            BigFloat[] betanB = new BigFloat[N];
            betanB[0] = -1.52763883147;
            betanB[1] = 1.048327004372e-1;
            betanB[2] = 2.669121419134012e-2;

            // f vector
            BigFloat[] fnB = new BigFloat[N];

            // jacobian matrix
            var mx = Matrix.Create<BigFloat>(N, N);

            double stepSize = (1 - (1 / (float)N)) / N;
            double[] zi = Range(1 / N, 1, stepSize).ToArray();
            double[] zj = new double[N];
            Array.Copy(zi, 1, zj, 0, N);

            for (int l = 0; l < count; l++)
            {
                for (int n = 0; n < N; n++)
                {
                    fnB[n] = f(zj[n], betanB);

                    BigFloat[] gradB = diffedRangeB(zj[n], betanB);

                    for (int m = 0; m < N; m++)
                    {
                        mx[n, m] = gradB[m];
                    }
                }

                var snB = mx.GetCholeskyDecomposition();

                //var sn = JnB.Decompose(false);
                var snA = snB.Solve(fnB);
                var betanp1 = betanB.Zip(snA, (d1, d2) => d1 - d2).ToArray();

                // Check for convergence
                BigFloat tolChk = snA.OneNorm();
                BigFloat tarVal = Math.Pow(10, -M);
                if (tolChk < tarVal)
                {
                    // compute alpha    
                    myalpha = 1 / g(1, betanp1);

                    Console.WriteLine(l + ": " + tolChk);
                    Console.WriteLine("\n Converged! with alpha = " + myalpha);
                    return false;
                }

                betanB = betanp1;
                Console.WriteLine(l + ": " + tolChk);// + "; " + (1 / gD((double)1, betanp1)));
            }

            return true;
        }


        public static IEnumerable<double> Range(double min, double max, double step)
        {
            double i;
            for (i = min; i <= max; i += step)
                yield return i;

            if (i != max + step) // added only because you want max to be returned as last item
                yield return max;
        }


        private BigFloat g(BigFloat z, BigFloat[] a)
        {
            BigFloat s = 0.0f;
            for (int i = a.Length; i >= 1; i--)
            {
                s = z * z * (s + a[i]);
            }
            return 1 + s;
        }

        private BigFloat f(BigFloat z, BigFloat[] a)
        {
            BigFloat one = new BigFloat(1.0);
            BigFloat alpha = one / g(one, a);
            BigFloat r = g(z, a) - alpha * g(g(z / alpha, a), a);

            return r;
        }


        BigFloat[] diffedRangeB(BigFloat z, BigFloat[] a)
        {
            Solver<BigFloat> bigFloatSolver = new Solver<BigFloat>();
            Func<BigFloat[], BigFloat> fbeta3 = x => f(z, a);
            //bigFloatSolver.TargetFunction = fbeta3;

            //var fDiffs = new FiniteDifferences(N, fbeta3);
            //return fDiffs2.Gradient(a);
            return a;
        }



        double compute_delta()
        {
            var maxIt = 16;
            var maxItJ = 10;
            var a1 = 1.0;
            var a2 = 0.0;
            var d1 = 3.2;
            Console.WriteLine(" i       d");
            for (int i = 2; i <= maxIt; i++)
            {
                var a = a1 + (a1 - a2) / d1;
                for (int j = 1; j <= maxItJ; j++)
                {
                    var x = 0.0;
                    var y = 0.0;
                    for (int k = 1; k <= 1 << i; k++)
                    {
                        y = 1.0 - 2.0 * y * x;
                        x = a - x * x;
                    }
                    a -= x / y;
                }
                var d = (a1 - a2) / (a - a1);
                Console.WriteLine("{0,2:d}    {1:f8}", i, d);
                d1 = d;
                a2 = a1;
                a1 = a;
            }
            return d1;
        }


    }


}
