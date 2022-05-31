import math

class WeightedPermutations:

    "Get all lists of genetic permutations for n people with k genes"

    def __init__(self, n:int, k:int):

        self.permutations = []
        base_array = (1., 1., [None]*n)

        self.__generateWeightedStrings(n, base_array, 0, k)

    @staticmethod
    def run(n, k=1):

        temp = WeightedPermutations(n, k)

        return temp.permutations

    def __generateWeightedStrings(self, n, arr, i, k):

        if i==n:
            self.permutations += [arr]
            return None

        # Get all arrays with j in the i-th position
        # Position 0 is incremented by weight
        for j in range(k+1):
            temp_arr = arr[2].copy()
            temp_arr[i] = j
            new_weight = arr[0] * math.comb(k, j)
            self.__generateWeightedStrings(n, (new_weight, new_weight, temp_arr), i+1, k)

        return 'Complete'