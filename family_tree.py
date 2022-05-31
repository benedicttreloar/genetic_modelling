import logging
logger = logging.getLogger(__name__)
logger.setLevel('INFO')

class FamilyTree:

    def __init__(self):

        ''' Set up persistent objects to populate data '''
        self._family_list = []      # Numbered list of family members
        self.__relationships = []   # List of parent/child pairs
        self.__known_statuses = []      # List of known statuses
        self.recessive_list = None  # List of individuals known to be recessive carriers
        self._simulation = None     # Set simulation object

    @property
    def family_list(self):
        return self._family_list

    @property
    def relationships(self):
        return self.__relationships

    def get_status_dict(self, all=False):
        if all:
            return {name: status for (name, status) in self.__known_statuses}
        else:
            raise NotImplementedError

    @property
    def independent_genomes(self):

        ''' Get the list of people in the tree with no independent genomes '''
        return [person for person in self._family_list if person not in
                [child for (parent, child) in self.__relationships]]

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
                    self.__relationships += [(parent, child)]
                elif self.__check_person(child):
                    logger.error('Parent {} does not exist - add to family tree'.format(parent))
                else:
                    logger.error('Child {} does not exist - add to family tree'.format(child))

    def add_status(self, person:str, status:int):

        '''Store information about genetic status, checking is that person exists'''
        assert(self.__check_person(person))
        assert(status in [0,1,2])
        self.__known_statuses += [(person, status)]
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