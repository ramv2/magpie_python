import numpy as np
class EqualSumCombinations:
    """
    Class to generate all combinations of non-negative integers that have
    equal sum.
    """
    def __init__(self, sum, size):
        """
        Constructor to initialize the variables.
        :param sum: Desired sum (must be greater than 0).
        :param size: Number of integers (must be greater than 1).
        """
        if sum <= 0:
            raise ValueError("Sum must be positive.")
        if size < 2:
            raise ValueError("Size must be greater than 1.")

        self.sum = sum
        self.size = size
        self.comb = np.zeros(size)
        self.comb[0] = sum

    def get_combinations(self):
        """
        Function to generate all the combinations of non-negative integers
        that have equal sum.
        :return output: A list of numpy arrays representing the different
        possible combinations.
        """
        output = []
        while True:
            output.append(self.comb.copy())
            self.increment_counter()
            if self.comb[0] == self.sum:
                break
        return output

    def increment_counter(self):
        """
        Function to increment a counter that iterates over all possible
        combinations. States at [0, ..., sum] and ends at [sum, ..., 0].
        :return:
        """
        for pos in xrange(self.size-2, -1, -1):
            if self.comb[pos] > 0:
                self.comb[pos] -= 1
                self.comb[pos + 1] += 1
                return
            else:
                self.comb[pos] = self.comb[pos + 1]
                self.comb[pos + 1] = 0

if __name__ == "__main__":
    x = EqualSumCombinations(5, 2)
    for i in x.get_combinations():
        print i