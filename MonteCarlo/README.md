# Monte Carlo Simulations

This sub-repo contains python code for monte carlo simulation methods. Ensure the following Python packages are installed before using.
NOTE: The version of the package is not particular (this was developed using the latest package versions).
- math
- pandas
- numpy
- scipy

## MonteCarloEngine.py Usage
The monte carlo engine is a framework for conducting monte carlo simulations on the portfolio level, using correlation matrices and the Cholesky Factorization technique.
Simulated price data is based on the following equation:

current price = previous price * e^(sigma * sqrt(1) * current transformed return)

Where:
sigma = standard devation of the actual time-series returns
current transformed return = dot product of the randomly sampled returns and the factorization matrix

mcEngine = MonteCarloEngine(priceData, numObservations)
simData = mcEngine.Simulate()
simData = mcEngine.Simulate(symbol, nSims)

NOTE: adding the parameters will return a dataframe where the columns are the sequence of simulation iterations.

With providing the parameters in the Simulate function such as:
symbol='AAPL'
nSims=100
One can create the infamous monte carlo simulation graph:

![AAPL Sim](https://github.com/tzabcoder/FinancialRiskManagement/assets/60833046/e6d11d5d-eabb-4f95-8d79-a48f4dfbc0c9)
