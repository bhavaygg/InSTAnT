import numpy as np
from scipy.sparse import csr_matrix

class PoissonBinomial():
    def __init__(self, probabilities):
        """Initialize the class and calculate the ``pmf``.
        :param probabilities: sequence of success probabilities :math:`p_i \\in
            [0, 1] \\forall i \\in [0, N]` for :math:`N` independent but not
            identically distributed Bernoulli random variables
        :type probabilities: numpy.array
        """
        self.success_probabilities = np.array(probabilities)
        self.number_trials = len(self.success_probabilities)
        self.pmf_list = self.pmf_all()

    def pmf_all(self):
        curr_n_pmf = np.zeros(self.number_trials + 1)
        curr_n_pmf[0] = 1. - self.success_probabilities[0]
        curr_n_pmf[1] = self.success_probabilities[0]

        for i in range(1,self.number_trials):   #i is actually (i+1)-th trial
            curr_n_pmf = self.pmf_ith_row_simple(curr_n_pmf, i+1)

        return curr_n_pmf

    def pmf(self, num_success):
        return self.pmf_list[num_success]

    def pmf_ith_row_simple(self, curr_n_pmf, curr_n):
        #first curr_n is 2
        next_n_pmf = np.zeros(curr_n_pmf.shape)
        next_n_pmf[0] = curr_n_pmf[0] * (1 - self.success_probabilities[curr_n-1])

        for k in range(1,curr_n+1):
            next_n_pmf[k] = curr_n_pmf[k-1] * self.success_probabilities[curr_n-1] + curr_n_pmf[k] * (1 - self.success_probabilities[curr_n-1])
        return next_n_pmf


    def pmf_ith_row(self, curr_n_pmf, curr_n):
        curr_matrix = self.transition_matrix(curr_n)
        return np.matmul(curr_matrix, curr_n_pmf)

    def pval(self, num_success):   #One sided p val
        return self.pmf_list[num_success:].sum()

    def transition_matrix(self, i):
        temp_diag = np.diag([1. - self.success_probabilities[i-1]] * (i+1) + [0] * (self.number_trials-i))
        temp_diag_1 = np.diag([self.success_probabilities[i-1]] * i + [0] * (self.number_trials-i), 1)
        return temp_diag + temp_diag_1.T






