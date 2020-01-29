using System;
using System.Collections.Generic;
using System.Linq;
using DiffSharp.AD.Forward;

using Accord.Math;
using Accord.Math.Differentiation;

namespace FeigenbaumAlphaCalc
{
    class Program
    {
        int N = 700;
        int M = 1000000;
        double myalpha = 0;
        int count = 10;

        static void Main(string[] args)
        {            
            Program P = new Program();
            P.compute_alpha();
        }

        bool compute_alpha()
        {
            double[] betan = new double[N];
            betan[0] = -1.52763883147;
            betan[1] = 1.048327004372e-1;
            betan[2] = 2.669121419134012e-2;

            // f vector
            double[] fn = new double[N];

            // jacobian matrix
            double[,] Jn = new double[N, N];

            double stepSize = (1 - (1 / (float)N)) / N;
            double[] zi = Range(1 / N, 1, stepSize).ToArray();
            double[] zj = new double[N];
            Array.Copy(zi, 1, zj, 0, N);

            for (int l = 0; l < count; l++)
            {
                for (int n = 0; n < N; n++)
                {
                    fn[n] = fD(zj[n], betan);
                    double[] grad = diffedRange(zj[n], betan);

                    for (int m = 0; m < N; m++)
                    {
                        Jn[n, m] = grad[m];
                    }
                }

                var sn = Jn.Decompose(false);
                var snA = sn.Solve(fn);
                var betanp1 = betan.Zip(snA, (d1, d2) => d1 - d2).ToArray();

                // Check for convergence
                double tolChk = snA.ToMatrix().Norm1();
                double tarVal = Math.Pow(10, -M);
                if (tolChk < tarVal)
                {
                    // compute alpha    
                    myalpha = 1 / gD((double)1, betanp1);

                    Console.WriteLine(l + ": " + tolChk);
                    Console.WriteLine("\n Converged! with alpha = " + myalpha);
                    return false;
                }

                betan = betanp1;
                Console.WriteLine(l + ": " + tolChk);
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

        /*
        DiffSharp.Interop.D gD(DiffSharp.Interop.D z, double[] a)
        {
            DiffSharp.Interop.D s = DiffSharp.Interop.D.Zero;
            
            for (int i = a.Length - 1; i >= 0; i--)
            {
                s = z * z * (s + a[i]);
            }
            return 1 + s;
        }


        private DiffSharp.Interop.D fD(DiffSharp.Interop.D z, double[] a)
        {
            DiffSharp.Interop.D one = DiffSharp.Interop.D.One;
            DiffSharp.Interop.D alpha = one / gD(one, a);
            DiffSharp.Interop.D r = gD(z, a) - alpha * (gD(gD(z / alpha, a), a));

            return r;
        }

        D gD(D z, double[] a)
        {
            D s = D.Zero;

            for (int i = a.Length - 1; i >= 0; i--)
            {
                s = z * z * (s + a[i]);
            }
            return 1 + s;
        }

        private D fD(D z, double[] a)
        {
            D one = D.One;
            D alpha = one / gD(one, a);
            D r = gD(z, a) - alpha * (gD(gD(z / alpha, a), a));

            return r;
        }
        */

        double gD(double z, double[] a)
        {
            double s = 0;
            for (int i = a.Length - 1; i >= 0; i--)
            {
                s = z * z * (s + a[i]);
            }
            return 1 + s;
        }


        private double fD(double z, double[] a)
        {
            double one = 1;
            double alpha = one / gD(one, a);
            double r = gD(z, a) - alpha * (gD(gD(z / alpha, a), a));
            return r;
        }

        private D gradfD(double z, double[] a)
        {
            D result;
            double estimatedError;
            D zD = D.NewD(z);

            var xvec = new MathNet.Numerics.LinearAlgebra.Double.DenseVector(a.Length);
            var yvec = new MathNet.Numerics.LinearAlgebra.Double.DenseVector(a.Length);
            var dydx = new MathNet.Numerics.LinearAlgebra.Double.DenseVector(a.Length);

            //Func<double, double> fbeta = x => fD(z, a);
            //Func<D, D[], D> fbeta2 = (x, y) => fD(zD, a);



            DiffSharp.Interop.D[] aD = new DiffSharp.Interop.D[100];
            D[] aDF = new D[100];

            // convert double arr to D[], two types
            // also put in to densearrays  
            for (int i = 0; i < a.Length; i++)
            {
                aD[i] = a[i];                  // interop type
                //aDF[i] = fbeta2(zD, aDF);      // forward type

            }

            // nonsense
            //result = fbeta2(zD, aDF);

            //result = fbeta.ForwardDerivative(a[0], out estimatedError);


            //var gfDel = fbeta2(zD, aDF);
            
            //var gf = DiffSharp.Interop.AD.Grad(x => fbeta(z, a), a);
            //var gf1 = DiffOps.grad(fbeta2(zD, aDF), aDF);

            //var iRes = gf.DynamicInvoke(aD);
            return D.One;
        }

        double[] diffedRange(double z, double[] a)
        {
            Func<double[], double> fbeta3 = x => fD(z, a);
            var fDiffs = new FiniteDifferences(N, fbeta3);
            double[] result = fDiffs.Gradient(a);
            return result;
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
