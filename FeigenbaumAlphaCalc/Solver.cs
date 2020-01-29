using System;
using Extreme.Mathematics.Generic;
using Extreme.Mathematics;


namespace FeigenbaumAlphaCalc
{
    /// <summary>
    /// Class that contains the generic Newton-Raphson algorithm.
    /// </summary>
    /// <typeparam name="T">The operand type.</typeparam>
    class Solver<T>
    {
        // Use a static variable to hold the arithmetic instance.
        // We get the instance from the global TypeAssociationRegistry:
        static IFieldOperations<T> ops =
            (IFieldOperations<T>)TypeAssociationRegistry.GetInstance(
                typeof(T), TypeAssociationRegistry.ArithmeticKey);

        // Member fields:
        Func<T, T> f, df;
        int maxIterations = 100;

        // The function to solve:
        public Func<T, T> TargetFunction
        {
            get { return f; }
            set { f = value; }
        }
        // The derivative of the function to solve.
        public Func<T, T> DerivativeOfTargetFunction
        {
            get { return df; }
            set { df = value; }
        }
        // The maximum number of iterations.
        public int MaxIterations
        {
            get { return maxIterations; }
            set { maxIterations = value; }
        }

        // The core algorithm.
        // Arithmetic operations are replaced by calls to
        // methods on the arithmetic object (ops).
        public T Solve(T initialGuess, T tolerance)
        {
            int iterations = 0;

            T x = initialGuess;
            T dx = ops.Zero;
            do
            {
                iterations++;
                // Compute the denominator of the correction term.
                T dfx = df(x);
                // Relational operators map to the Compare method.
                // We also use the value of zero for the operand type.
                // if (dfx == 0)
                if (ops.Compare(dfx, ops.Zero) == 0)
                {
                    // Change value by 2x tolerance.
                    // When multiplying by a power of two, it's more efficient 
                    // to use the ScaleByPowerOfTwo method.
                    dx = ops.ScaleByPowerOfTwo(tolerance, 1);
                }
                else
                {
                    // dx = f(x) / df(x)
                    dx = ops.Divide(f(x), dfx);
                }
                // x -= dx;
                x = ops.Subtract(x, dx);

                // if |dx|^2<tolerance
                // Convergence is quadratic (in most cases), so we should be good here:
                if (ops.Compare(ops.Multiply(dx, dx), tolerance) < 0)
                    return x;
            }
            while (iterations < MaxIterations);
            return x;
        }
    }
}

