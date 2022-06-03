# Carrier Probability for Recessive Genetic Disorders

Estimate carrier probabilities within family trees of known or unknown recessive genetic characteristics

This code allows the user to:

- define the structure of a family tree
- specify known genetic statuses for various family members
- simulate all potential transmission chains for all valid starting conditions
- given the prior population likelihood of a given genetic status, estimate the carrier status probability given the structure of the family tree and the known constraints

For example, this allows investigation of the prevalence and frequency of blue-eyed genes within a family structure, or target carrier screening for severe recessive genetic conditions such as Mucopolysaccharidosis 

### Example Notebook

For an example genetic tree, use:

_genetic_analysis_public.ipynb_

### Structure

##### FamilyTree Object

This allows the user to define the structure of a given family tree.

The principal method is _add_member_, which takes a person's name plus optionally the identity of their parents and their known carrier status.
The parents variable can be a string or a list and requires the parents to already be defined. 
The status variable can currently take values between 0 and 2, assuming a single gene recessive status.
This method calls a further two methods, _add_relationship_ which allows the user to define parental relationships, and _add_status_ which lets people define their genetic statuses. 

The object then has a number of properties which can be called by any simulation module to define the structure of the tree and return known genetic statuses.

##### SimulatePrevalence Object

This takes a defined family tree and simulates all genetic transmission chains to provide an estimate of carrier probability. 
It also provides an estimate of how much effect a test would have on the probabilities in the family tree.

Most methods are hidden, with the simulation run as part of the init call. 
The relevant method is _individual_probabilities_, which return the carrier status probability for each person.

##### WeightedPermutations

This provides a list containing all possible transmission permutations, weighted by the frequency with which they would occur.

### Next Steps

Potential next steps:

- The tree structure should be optimised to ensure it stops as early as possible in the modelling of the transmission chain to increase speed
- Adding metrics showing how much effect a test for each person will have on the conditional carrier probabilities - this provides a metric for the value of testing specific members
- Switching to a Bayesian modelling framework with constraints from a simulation framework, which should allow for much faster calculation and enable moving towards more complex, multi-site genetic conditions such as ADHD 