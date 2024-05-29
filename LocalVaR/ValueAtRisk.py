# File Imports
import numpy as np

"""
* VaR - Value at Risk
*
* Descritption:
* This class calculates the local value at risk (VaR) for any given portfolio
* and its corresponding weights in percentage terms.
* The VaR class calculates the total VaR, individual VaR, marginal VaR,
* and component VaR.
"""
class VaR:
    _stddevMatrix = []
    _corrMatrix = []
    _covMatrix = []
    
    # Confidence interval variables
    # NOTE: each index of the confidence intervals associates with
    #       the corresponding z-score
    _validConfInts = [
        0.90,  # 90% confidence interval
        0.91,  # 91% confidence interval
        0.92,  # 92% confidence interval
        0.93,  # 93% confidence interval
        0.94,  # 94% confidence interval
        0.95,  # 95% confidence interval
        0.96,  # 96% confidence interval
        0.97,  # 97% confidence interval
        0.98,  # 98% confidence interval
        0.99,  # 99% confidence interval
        0.995  # 99.5% confidence interval
    ]
    _zscores = [
        1.2816,
        1.3408,
        1.4051,
        1.4758,
        1.5548,
        1.6449,
        1.7507,
        1.8808,
        2.0537,
        2.3263,
        2.5758
    ]

    """
    * __init__()
    *
    * Initializes the VaR instance by verifying the construction params
    * and constructing the required matrices.
    * @param[in] returnMatrix(list)   - matrix of daily returns [[], [], ...]
    * @param[in] positionMatrix(list) - list of portfolio weights in decimal percentages
    """
    def __init__(self, returnMatrix, positionMatrix):
        # Check the type of the returnMatrix and positionMatrix
        if returnMatrix is None or type(returnMatrix) != list:
            print('returnMatrix must be a list of returns...')
        
        elif positionMatrix is None or type(positionMatrix) != list:
            print('positionMatrix must be a list of postitions...')

        else:
            # Set the local return and position matrix
            self._returnMatrix = returnMatrix
            self._positionMatrix = positionMatrix

            self._ConstructVarMatrix()

    """
    * _ConstructVarMatrix()
    *
    * Calculate the required matrices for the VaR calculations.
    """
    def _ConstructVarMatrix(self):
        # Calculate the standard deviations
        for returns in self._returnMatrix:
            self._stddevMatrix.append(np.std(returns))

        # Calculate the correlation matrix
        self._covMatrix = np.cov(self._returnMatrix)

        # Calculate the covariance matrix
        self._corrMatrix = np.corrcoef(self._returnMatrix)

        # Calculate the coeffecient matrix
        self._coeffMatrix = np.dot(self._covMatrix, self._positionMatrix)
        self._posCoeffMatrix = list(np.multiply(self._positionMatrix, self._coeffMatrix))
    
    """
    * _GetZscore()
    *
    * Get the Z-score associated with a valid confidence interval.
    * @param[in] conf(float) - a valid confidence interval
    """
    def _GetZscore(self, conf):
        # Return the z-score at the corresponding confidence interval
        return self._zscores[self._validConfInts.index(conf)]

    """
    * TotalVaR()
    *
    * Calculates the total value at risk for the portfolio.
    * @param[in] confidence(float) - requested confidence interval for VaR
    """
    def TotalVaR(self, confidence=0.99):
        t_var = None

        # Verify confidence interval
        if confidence in self._validConfInts:
            _zScore = self._GetZscore(confidence)

            _total = sum(self._posCoeffMatrix) # position * coefficientMatrix
            _risk = np.sqrt(_total)

            t_var = _risk * _zScore

        return t_var
    
    """
    * IndividualVaR()
    *
    * Calculates the value at risk for each position in the portfolio.
    * @param[in] confidence(float) - requested confidence interval for VaR
    """
    def IndividualVaR(self, confidence=0.99):
        i_var = []

        # Verify the confidence interval
        if confidence in self._validConfInts:
            _zScore = self._GetZscore(confidence)

            # VaR for each position
            for i in range(len(self._positionMatrix)):
                i_var.append(_zScore * self._stddevMatrix[i] * self._positionMatrix[i])

        return i_var

    """
    * MarginalVaR()
    *
    * Calculates the marginal value at risk for each position in the portfolio.
    * @param[in] confidence(float) - requested confidence interval for VaR
    """
    def MarginalVaR(self, confidence=0.99):
        m_var = []

        # Verify confidence interval
        if confidence in self._validConfInts:
            _tVar = self.TotalVaR(confidence)

            _total = sum(self._posCoeffMatrix)
            _betaMatrix = [c/_total for c in self._coeffMatrix] # calculate beta for each pos
            m_var = [beta * _tVar for beta in _betaMatrix]      # mVaR = tVaR * beta

        return m_var

    """
    * ComponentVaR()
    *
    * Calculates the component value at risk for each position in the portfolio.
    * @param[in] confidence(float) - requested confidence interval for VaR
    """
    def ComponentVaR(self, confidence=0.99):
        c_var = []

        # Verify confidence interval
        if confidence in self._validConfInts:
            m_var = self.MarginalVaR(confidence)

            c_var = list(np.multiply(m_var, self._positionMatrix)) # cVaR = mVaR * pos

        return c_var
