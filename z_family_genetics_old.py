import logging
from weighted_permutations import WeightedPermutations

logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')

class FamilyGenetics:

    ''' Note: this class is unoptimised and contains both simulation and tree definition '''

    def __init__(self):
        # Numbered list of family members
        self._family_list = []

        # List of parent/child pairs
        self._family_tree = []

        # List of tested statuses
        self._status_list = []
        self.recessive_list = None

        # Set simulation object
        self._simulation = None

    @property
    def independent_genomes(self):

        ''' Get the list of people in the tree with no independent genomes '''

        return [person for person in self._family_list if person not in
                [child for (parent, child) in self._family_tree]]

    def add_member(self, person:str, parents:list=None, status:int=None):

        ''' Add a person to the family tree and add any included parental and status information '''

        # Check family member not yet added
        if not self.__check_person(person):
            self._family_list += [person]
        else:
            raise ValueError('Person {} already exists - disambiguate name'.format(person))

        self.add_relationship(person, parents)

        if not status is None:
            self.add_status(person, status)

    def add_relationship(self, child:str, parents:list):

        '''Store information about a parent/child relationship'''

        if parents is not None:
            if type(parents) == str:
                parents = [parents]
            for parent in parents:
                if self.__check_person(child) and self.__check_person(parent):
                    self._family_tree += [(parent, child)]
                elif self.__check_person(child):
                    logger.error('Parent {} does not exist - add to family tree'.format(parent))
                else:
                    logger.error('Child {} does not exist - add to family tree'.format(child))

    def add_status(self, person:str, status:int):

        '''Store information about genetic status, checking is that person exists'''

        assert(self.__check_person(person))
        assert(status in [0,1,2])
        self._status_list += [(person, status)]
        if status == 2.:
            self.recessive_list = self.__update_recessive_list(person)

    def __update_recessive_list(self, recessive_list:list):

        ''' Update the list of recessive carriers and validate status lists '''

        if recessive_list is None:
            return self.recessive_list
        elif type(recessive_list) == str:
            recessive_list = [recessive_list]

        # If there is already a recessive list, retrieve and append
        if self.recessive_list is not None:
            recessive_list = list(set(self.recessive_list + recessive_list))

        # Check each recessive carrier exists
        for person in recessive_list:
            assert(self.__check_person(person))

        return recessive_list

    def specify_recessive_list(self, recessive_list, genomes=None, name_list=None):

        ''' Specify a list of symptomatic carriers of a recessive gene and update status list '''

        # Store and update recessive list
        self.recessive_list = self.__update_recessive_list(recessive_list)

        # If no genetic information is provided use full simulation set
        if genomes is None:
            # Use simulation data as the input
            genomes = self._simulation
            name_list = self._family_list

        if recessive_list is None:
            return 0

        # Define new name list for speed
        new_name_list = recessive_list + [x for x in name_list if x not in recessive_list]

        # Enforce all people not in dominant list are not
        for person in name_list:
            if person in recessive_list:
                # Enforce all people in recessive list aren't showing the gene
                genomes = self._enforce_status((person, 2), genomes, name_list)
            else:
                # Enforce all people missing from recessive list aren't showing the gene
                genomes = self._enforce_status((person, 2), genomes, name_list, include=False)

        return genomes

    def _enforce_status(self, person_status:'(person, status)',
                        genome_in=None, names_in=None,
                        include:bool=True, keep=False):

        ''' Enforce status conditions on a set of genetic conditions '''

        (person, status) = person_status
        if genome_in is None:
            genome = self._simulation
            names_in = self._family_list
        else:
            genome = genome_in

        if person not in names_in:
            return genome

        if genome is not None:
            sum_wgt_before = sum([new_wgt for (prior_wgt, new_wgt, chain) in genome])
            if keep:
                # This option preserves all data; may pose issues with storage size and speed
                if include:
                    genome_out = [(prior_wgt, new_wgt, chain)
                                  if chain[names_in.index(person)] == status
                                  else (prior_wgt, 0., chain)
                                  for (prior_wgt, new_wgt, chain) in genome]
                else:
                    genome_out = [(prior_wgt, new_wgt, chain) if chain[names_in.index(person)] != status
                                  else (prior_wgt, 0., chain)
                                  for (prior_wgt, new_wgt, chain) in genome]
            else:
                # This option discards data which is wrong
                if include:
                    genome_out = [(prior_wgt, new_wgt, chain)
                                  for (prior_wgt, new_wgt, chain) in genome
                                  if chain[names_in.index(person)] == status]
                else:
                    genome_out = [(prior_wgt, new_wgt, chain)
                                  for (prior_wgt, new_wgt, chain) in genome
                                  if chain[names_in.index(person)] != status]

            sum_wgt_after = sum([new_wgt for (prior_wgt, new_wgt, chain) in genome_out])
            if sum_wgt_before > 0:
                logger.debug('Condition {} - {} is {}: {} of probability space dropped'.format(
                    person, status, include, 1 - sum_wgt_after/sum_wgt_before))

            return genome_out
        else:
            return None

    def __check_person(self, person):

        # Check person exists in tree
        check = (len([1 for x in self._family_list if x==person])==1)

        return check == 1

    def simulate(self, gene_frequency=1/300):

        #logger.info('Starting simulation')
        print('Starting simulation')
        # Get people with independent genomes and generate a set of genomes
        key_people = self.independent_genomes
        raw_initial_genomes = WeightedPermutations.run(len(key_people), 2)
        logger.info('{} key people; initially {} independent initial genomes'.format(
            len(key_people), len(raw_initial_genomes)))
        raw_initial_genomes = self._check_genomes(raw_initial_genomes, key_people)

        # Adjust initial genomes by prior distribution of population prevalence
        def weight_prevalence(transmission, gene_frequency):

            (weight, junk, chain) = transmission
            weight = weight * (gene_frequency ** sum(chain)) *\
                     ((1-gene_frequency) ** (len(chain) - sum(chain)))

            return (weight, weight, chain)

        initial_genomes = [weight_prevalence(transmission, gene_frequency)
                           for transmission in raw_initial_genomes]

        logger.info('{} key people; after adjustment {} independent initial genomes'.format(
            len(key_people), len(initial_genomes)))

        # Get potential transmission trees
        transmission_chains = WeightedPermutations.run(len(self._family_tree), 1)

        logger.info('{} potential transmission chains'.format(
            len(transmission_chains)))

        # Calculate actual genetic outcomes
        def simulate_genetics(input_statuses:dict, transmission_list:list):
            status = {name: 0 for name in self._family_list}
            for person in input_statuses.keys():
                status[person] = input_statuses[person]

            for pos in range(len(transmission_list)):
                (parent, child) = self._family_tree[pos]
                if status[parent] == 2:
                    status[child] = 1 + int(status[child])
                elif status[parent] * transmission_list[pos] == 1:
                    status[child] = int(status[child]) + status[parent] * transmission_list[pos]
                elif status[parent] == 0:
                    status[child] = int(status[child])

            return [status[name] for name in self._family_list]

        # Simulate genomes
        logger.info('Starting simulation - {} genomes to create'.format(
            len(initial_genomes) * len(transmission_chains)))
        genomes = [] # This will be: (prior_probability, new_probability, genome)
        for (orig_swgt, seed_wgt, seed) in initial_genomes:
            for (orig_cwgt, chain_wgt, chain) in transmission_chains:
                genome_wgt = seed_wgt * chain_wgt
                seed_status = dict(zip(key_people, seed))
                genome = simulate_genetics(seed_status, chain)

                genomes += [(genome_wgt, genome_wgt, genome)]

        logger.info('Simulation complete - checking {} genome sets'.format(
            len(genomes)))
        self._simulation = self._check_genomes(genomes, self._family_list)

        logger.info('Checks complete - {} genome sets remaining'.format(
            len(self._simulation)))

        return self._simulation

    def _check_genomes(self, genomes, name_list, keep=False):

        # Check a list of genomes against the set of currently stored conditions
        # Remove excluded outcomes
        genomes = [(prior_wgt, new_wgt, chain)
                   for (prior_wgt, new_wgt, chain) in genomes
                   if prior_wgt != 0.]
        for status in self._status_list:
            genomes = self._enforce_status(status, genomes, name_list, include=True, keep=keep)
        genomes = self.specify_recessive_list(self.recessive_list, genomes, name_list)

        return genomes

    def get_probabilities(self, status):

        probabilities = {}
        for person in self._family_list:
            pos = self._family_list.index(person)
            try:
                prob = sum([new for (old, new, chain) in self._simulation if chain[pos] == status])\
                       / sum([new for (old, new, chain) in self._simulation])
            except:
                prob=False
            probabilities[person] = prob
            print('{}: {}'.format(person, prob))

        return probabilities