# File Imports
import math
import pandas as pd
import numpy as np
import scipy.linalg
from scipy import special

"""
* Monte Carlo Simulation
*
* Description:
* This class is an engine for the monte carlo simulation.
* Given initial conditions and a group of timeseries data, the class will generate (simulate)
* price data for a given set of assets. The simulation technique is based on the cholesky
* decomposition and random return transformation method.
"""
class MonteCarloEngine:
    # Class Variables
    _simulatedPriceData = None

    """
    * __init__() : public
    *
    * Verify input parameter types and set the engine parameters.
    * Constructor than creates the actual returns for the price data for
    * the provided price matrix.
    * @param[in] priceMatrix(pd.DataFrame) - dataframe of timeseries prices
    * @param[in] observations(int)         - number of periods into the future to simulate the price data
    """
    def __init__(self, priceMatrix, observations=1000):
        if type(priceMatrix) != pd.DataFrame and type(observations) != int:
            print("The price matrix must be a dataframe")

        else:
            # Set the monte carlo parameters
            self._observations = observations
            self._priceMatrix = priceMatrix

            self._CreateReturnMatrix() # Calculate the returns from the actual price data

    """
    * _CreateReturnMatrix() : private
    *
    * Calculates the returns for each column in the priceMatrix. (Calculates the return
    * for each price timeseries). The returns are added to a seperate return matrix and
    * their standard deviations are calculated.
    *
    * NOTE: The column names correspond to the timeseried (priceMatrix) column names.
    """
    def _CreateReturnMatrix(self):
        self._returnMatrix = pd.DataFrame()
        self._stddevMatrix = []

        for col, data in self._priceMatrix.items():
            # Calculate returns and stddev
            self._returnMatrix[col] = data.pct_change()
            self._stddevMatrix.append(self._returnMatrix[col].std())

    """
    * _CreateRandomReturns() : private
    *
    * For each provided timeseries, create a list of randomly generated returns. The random
    * returns are calculated by randomly sampling the inverse of the standard normal distribution.
    * After returns are created, it generates a random return matrix.
    *
    * NOTE: The column names correspond to the timeseried (priceMatrix) column names.
    """
    def _CreateRandomReturns(self):
        _randomReturnMatrix = {}

        for col in self._returnMatrix.columns:
            # Randomly sample the inverse of the standard normal distribution
            _tempRandomReturns = [special.ndtri(np.random.rand()) for i in range(self._observations)]
            _randomReturnMatrix[col] = _tempRandomReturns

        # Create the random matrix from the generated data
        self._randomMatrix = pd.DataFrame.from_dict(_randomReturnMatrix)

    """
    * _CreateCorrelationMatrix() : private
    *
    * Creates the correlation matrix for the random return matrix.
    *
    * * NOTE: The column/row names correspond to the timeseried (priceMatrix) column names.
    """
    def _CreateCorrelationMatrix(self):
        self._correlationMatrix = self._randomMatrix.corr()

    """
    * _CholeskyDecomposition() : private
    *
    * A = LL^T
    * Calculates the cholesky foactorization for the correlation matrix. The correlation
    * matrix is factored into L and U(L^T). The factorization matrix is the transposed L.
    """
    def _CholeskyDecomposition(self):
        A = np.array(self._correlationMatrix)
        L = scipy.linalg.cholesky(A, lower=True)
        U = scipy.linalg.cholesky(A, lower=False)

        self._factorizationMatrix = U

    """
    * _TransformRandomReturns() : private
    *
    * Transforms the random returns by multiplying (dot product) the random returns and the
    * factorization matrix from the cholesky factorization.
    * The function creates a trandformed returns matrix.
    *
    * NOTE: The column names correspond to the timeseried (priceMatrix) column names.
    """
    def _TransformRandomReturns(self):
        # Multiply the random returns buy the factorization matrix
        _transformed = np.dot(self._randomMatrix, self._factorizationMatrix)

        # Convert transformed returns into a dataframe
        self._transformedReturns = pd.DataFrame(_transformed, columns=self._returnMatrix.columns)

    """
    * _SimulatePriceData() : private
    *
    * Simulates the price data for each asset initial asset. After calculating the
    * simulated price, the data is compiled into a data matrix.
    *
    * simulatedPrice(i) = price(i-1) * e^(sigma * sqrt(1) * transformedReturn(i))
    *
    * NOTE: The column names correspond to the timeseried (priceMatrix) column names.
    """
    def _SimulatePriceData(self):
        _simulatedPrices = {}

        stddevIdx = 0 # index for tracking correct standard deviation
        for col, transRet in self._transformedReturns.items():
            # obtain initial price and standard deviation
            _initialPrice = self._priceMatrix[col].iloc[-1]
            _stddev = self._stddevMatrix[stddevIdx]

            _prices = []
            for i in range(len(transRet)):
                if i == 0:
                    _prices.append( _initialPrice * math.exp(_stddev * np.sqrt(1) * transRet.iloc[i]) )

                else:
                    _prices.append( _prices[i-1] * math.exp(_stddev * np.sqrt(1) * transRet.iloc[i]) )

            _simulatedPrices[col] = _prices
            stddevIdx += 1

        self._simulatedPriceData = pd.DataFrame.from_dict(_simulatedPrices)

    """
    * Simulate() : public
    *
    * Runs the monte carlo simulation for the input data.
    * @return _simulatedPriceData(pd.DataFrame)
    """
    def Simulate(self, symbol=None, nSims=1):
        self._CreateRandomReturns()     # Create random returns for each asset in the price matrix
        self._CreateCorrelationMatrix() # Calculate the correlations from the random returns
        self._CholeskyDecomposition()   # Decompose the correlation matrix through the Cholesky Factorization
        self._TransformRandomReturns()  # Transform returns by multiplying the random returns by the factorized matrix
        self._SimulatePriceData()       # Simulate the price data

        return self._simulatedPriceData