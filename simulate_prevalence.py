import pandas as pd
import logging
from weighted_permutations import WeightedPermutations
from family_tree import FamilyTree

logger = logging.getLogger(__name__)
logger.setLevel('INFO')


class SimulatePrevalence:

    def __init__(self, family_tree:FamilyTree,
                 recessives_are_known:bool=True, recessive_prevalence:float=None):

        ''' Set up persistent objects to populate data '''
        self._recessives_are_known=recessives_are_known
        self._recessive_prevalence=recessive_prevalence
        self.optimise_tree(family_tree)

        # Define simulation object
        self._simulation = None
        self.simulate(recessive_prevalence)# Set simulation object

    def optimise_tree(self, family_tree:FamilyTree):

        self._original_tree = family_tree
        self.__family_list = family_tree.family_list
        self.__relationships = family_tree.relationships
        self.__independent_genomes = family_tree.independent_genomes
        if family_tree.recessive_list is None:
            self.__recessive_list = None
        else:
            self.__recessive_list = set(family_tree.recessive_list)
        self._status_dict = family_tree.get_status_dict(all)

    def _get_possible_initial_genomes(self, gene_frequency):

        ''' Generate the set of initial independent genomes '''

        # Get the list of people with independent genomes
        key_people = self.__independent_genomes
        logger.debug('Key people are: {}'.format(', '.join(key_people)))

        # Generate all initial genomes
        raw_chains = WeightedPermutations.run(len(key_people), 2)
        logger.info('{} key people; initially {} independent unvalidated initial genomes'.format(
            len(key_people), len(raw_chains)))

        # Enforce conditions on initial genomes
        raw_dicts = [(weight1, weight2, dict(zip(key_people, chain)))
                     for (weight1, weight2, chain) in raw_chains]
        checked_initial_dicts = self._check_list_of_genomes(raw_dicts)

        # Adjust initial genomes by prior distribution of population prevalence
        def weight_prevalence(transmission, gene_frequency):

            (weight, junk, gene_dict) = transmission
            chain = list(gene_dict.values())
            weight = weight * (gene_frequency ** sum(chain)) *\
                     ((1-gene_frequency) ** (2 * len(chain) - sum(chain)))

            return (weight, weight, chain)

        weighted_genomes = [weight_prevalence(transmission, gene_frequency)
                           for transmission in checked_initial_dicts]

        def convert_to_dicts(weighted_genomes, key_people):

            def convert(potential_initial, key_people):

                (weight, weight, chain) = potential_initial

                return (weight, weight, dict(zip(key_people, chain)))

            return [convert(x, key_people) for x in weighted_genomes]

        logger.info('{} key people; after validation, {} independent initial genomes remain'.format(
            len(key_people), len(weighted_genomes)))

        return convert_to_dicts(weighted_genomes, key_people)

    def _get_transmission_chains(self, relationships=None):

        if relationships is None:
            relationships = self.__relationships

        logger.debug('Transmission chain simulation starting')
        transmissions = WeightedPermutations.run(len(relationships), 1)
        logger.debug('Transmission chain simulation finished')
        logger.info('{} potential transmission chains'.format(len(transmissions)))

        return transmissions

    def _check_list_of_genomes(self, genomes):

        temp = [self._check_genome(genome) for genome in genomes]

        return [x for x in temp if x is not None]

    def _check_genome(self, genome):

        (prior_wgt, new_wgt, genome_dict) = genome

        for person in genome_dict.keys():
            if not self._check_status_is_valid(genome_dict[person], person):
                return None
            if not self._check_recessive_list(genome_dict[person], person):
                return None

        return genome

    def _check_status_is_valid(self, status, child):

        ''' Check returned genetic status given chain is valid '''
        try:
            if child in self._status_dict.keys():
                possible_statuses = self._status_dict[child]
                if type(possible_statuses) is list:
                    return status in possible_statuses
                else:
                    return status == possible_statuses
            else:
                return True
        except:
            print(self._status_dict)
            print(child)
            raise ValueError

    def _check_recessive_list(self, status, child):

        ''' Check known recessive child in status list '''
        if self._recessives_are_known and (self.__recessive_list is not None):
            if status == 2:
                return (child in self.__recessive_list)
            else:
                return (child not in self.__recessive_list)
        else:
            return True

    def _create_genome_sets(self, initial_genomes, transmission_chains):

        ''' Create probability weighted sets of genomes

        Iterate through pairs of initial conditions and transmission chains '''

        genomes = []
        logger.info('Starting simulation - {} potential genome sets to create'.format(
            len(initial_genomes) * len(transmission_chains)))
        for (orig_swgt, seed_wgt, seed_status) in initial_genomes:
            for (orig_cwgt, chain_wgt, chain) in transmission_chains:
                genome_wgt = seed_wgt * chain_wgt
                genome = self._generate_single_chain(seed_status, chain)
                if genome is not None:
                    genomes += [(genome_wgt, genome_wgt, genome)]

        logger.info('Simulation complete - {} valid genome sets returned'.format(
            len(genomes)))

        return genomes

    def _generate_single_chain(self, input_statuses:dict, transmission_chain:list):

        ''' Generate a full set of genetic material given input genetics and a transmission list

        Returns: a tuple with (weight, weight, [genetic_status])
        '''

        # Create status dictionary
        status = {name: 0 for name in self.__family_list}

        # Enforce initial conditions
        for person in input_statuses.keys():
            status[person] = input_statuses[person]

        # Create indicator to check if child genome is complete
        child_list = [child for (parent, child) in self.__relationships]
        def check_pos(child_list, pos):
            try:
                return (child_list[pos] not in child_list[pos+1:])
            except:
                return False
        child_complete = [check_pos(child_list, pos) for pos in range(len(child_list))]

        # Generate transmission chain
        for pos in range(len(transmission_chain)):

            # Get implied genetic status
            (parent, child) = self.__relationships[pos]
            if status[parent] == 2:
                status[child] = 1 + int(status[child])
            elif status[parent] * transmission_chain[pos] == 1:
                status[child] = int(status[child]) + status[parent] * transmission_chain[pos]
            elif status[parent] == 0:
                status[child] = int(status[child])

            # Check status is valid
            if child_complete[pos]:
                if not self._check_status_is_valid(status[child], child):
                    return None
                if not self._check_recessive_list(status[child], child):
                    return None

        # Return an ordered list of genetic statuses
        return [status[name] for name in self.__family_list]

    def simulate(self, recessive_prevalence=None):

        #logger.info('Starting simulation')
        print('Starting simulation')

        # Get initial genomes to seed with
        if recessive_prevalence is None:
            recessive_prevalence = self._recessive_prevalence
        assert(recessive_prevalence is not None)
        gene_frequency = recessive_prevalence ** 0.5
        initial_genomes = self._get_possible_initial_genomes(gene_frequency)

        # Get potential transmission trees
        transmission_chains = self._get_transmission_chains(self.__relationships)

        # Simulate genomes
        self._simulation = self._create_genome_sets(initial_genomes, transmission_chains)

        return self._simulation

    def individual_probabilities(self, statuses=[0,1,2], genomes=None):

        if genomes is None:
            genomes = self._simulation

        if type(statuses) in (int, float):
            statuses = [statuses]

        probabilities = {}
        for status in statuses:
            # Calculate probability of each status for the list of people
            def prob(status, person):
                pos = self.__family_list.index(person)
                try:
                    prob = sum([new for (old, new, chain) in genomes if chain[pos] == status])\
                           / sum([new for (old, new, chain) in genomes])
                except:
                    prob=False
                return prob
            probabilities[status] = [prob(status, person) for person in self.__family_list]

        return pd.DataFrame(probabilities, index=self.__family_list)