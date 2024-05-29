# File imports
import numpy as np

"""
* VaR - Value At Risk
*
* Description:
* This class calculates the full valuation VaR for any given portfolio
* in dollar terms.
* This class uses the historical simulation method to estimate VaR.
"""
class VaR:
    _returnMatrix = []
    _portfolioMatrix = []

    """
    * __init__()
    *
    * Initializes the VaR instance by creating the return and portfolio matrices.
    * The initialization then simulates the portfolio over the past data period.
    * @param[in] dataMatrix(pd.DataFrame) - portfolio Adj Close price data
    * @param[in] positionMatrix(list)     - list of shares for each position
    *
    * NOTE: dataMatrix column positions correspond to the positionMatrix index
    *       dataMatrix[col1] => positionMatrix[pos1]
    *       ...
    *       dataMatrix[colN] => positionMatrix[posN]
    """
    def __init__(self, dataMatrix: pd.DataFrame, positionMatrix: list) -> None:
        idx = 0 # index for positions

        for column, data in dataMatrix.items():
            # calculate returns and total positions
            self._returnMatrix.append(data.pct_change().dropna())
            self._portfolioMatrix.append(positionMatrix[idx] * data[::-1][0]) # shares * recent price

            idx += 1

        self._SimulatePortfolio()

    """
    * _SimulatePortfolio()
    *
    * Simulates the portfolio over the historical period using current positions.
    * Provides a daily p&l over the historical period in dollar terms. The daily p&l is
    * calculated through taking the dot product of the daily returns and the current
    * positions.
    """
    def _SimulatePortfolio(self) -> None:
        self._simulationVector = []

        for i in range(len(self._returnMatrix[0])):
            _ret = [r[i] for r in self._returnMatrix]
            self._simulationVector.append(np.dot(_ret, self._portfolioMatrix))

        # reverse order to match dataMatrix
        self._simulationVector.reverse()

    """
    * FullValuationVaR()
    *
    * Uses the historical method to estimate the full valuation VaR for the given
    * portfolio, and the given confidence interval.
    * @param[in] confidence(float) - confidence interval to calculate the VaR
    * @return fullVaR(float)
    """
    def FullValuationVaR(self, confidence: float=0.99) -> float:
        fullVaR = None

        # Validate the confidence interval
        if confidence > 0 and confidence <= 1:
            # VaR(Index) = N * (1 - CI)
            # N = number of samples
            # CI = confidence interval
            nSamples = len(self._simulationVector)
            cIdx = int(round((nSamples*(1-confidence)), 0))

            _sortedSimulation = sorted(self._simulationVector)

            fullVaR = _sortedSimulation[cIdx-1]

        return fullVaR
