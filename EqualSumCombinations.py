import numpy as np
import cProfile
class EqualSumCombinations:
    """
    Class to generate all combinations of non-negative integers that have
    equal sum.
    """
    def __init__(self, sum_, size):
        """
        Constructor to initialize the variables.
        :param sum_: Desired sum (must be greater than 0).
        :param size: Number of integers (must be greater than 1).
        """
        if sum_ <= 0:
            raise ValueError("Sum must be positive.")
        if size < 2:
            raise ValueError("Size must be greater than 1.")

        # self.total = sum_
        # self.start = [0]*size
        # self.size_ = size
        # self.flag = False

        self.dp = np.zeros((sum_+1, size+1), dtype=object)
        self.combs = self.get_combinations(sum_, size)

    # def __iter__(self):
    #     self.start[0] = self.total
    #     self.flag = True
    #     return self
    #
    # def next(self):
    #     if self.flag:
    #         self.flag = False
    #         return self.start
    #     for pos in xrange(self.size_-2, -1, -1):
    #         if self.start[pos] > 0:
    #             self.start[pos] -= 1
    #             self.start[pos + 1] += 1
    #             break
    #         else:
    #             self.start[pos] = self.start[pos + 1]
    #             self.start[pos + 1] = 0
    #
    #     if self.start[0] == self.total:
    #         raise StopIteration
    #     return self.start

    def get_combinations(self, sum, n):
        """
        A recursive function to generate the list of all non-negative integer
        combinations
        of a given size that have a given sum.
        :param sum: Desired sum.
        :param n: Desired size.
        :return: A list containing the lists of combinations.
        """

        # Check if the list of lists has already been generated. If yes,
        # return the value.
        if self.dp[sum][n]:
            return self.dp[sum][n]

        # Initialize list and check special cases first.
        tmp_list = []
        if n == 1:
            tmp_list = [[sum]]
        elif sum == 0:
            tmp_list = [[0] * n]
        else:
            for i in xrange(sum, -1, -1):
                for l in self.get_combinations(sum - i, n - 1):
                    tmp_list.append([i]+l)

        # Store value into the numpy array before returning it.
        self.dp[sum][n] = tmp_list
        return tmp_list

if __name__ == "__main__":
    x = EqualSumCombinations(23, 7)
    cProfile.run('x.get_combinations(23, 7)')
    # l = x.get_combinations(23, 7)
    #
    # for i in xrange(7,24):
    #     print len(x.dp[i][6])
        # print len(x.get_combinations(i, 6))
    # for i,j in enumerate(EqualSumCombinations(5, 7)):
    #     print i,j