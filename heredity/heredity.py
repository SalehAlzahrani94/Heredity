import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """

    prob_of_gene=1
    for each_person in people: # go over each itme in people 
        #num means number
        if each_person in one_gene :# if the is one 
             gen_num=1 # save one gene
        elif each_person in two_genes:# if there is two 
             gen_num=2# do this part and save two  gene
        else :# else - there none
            gen_num=0# save two  gene 

        trait=True if each_person in have_trait else False
        gen_num_prob=PROBS['gene' ] [gen_num ]
        trait_prob=PROBS[ 'trait' ][ gen_num ][ trait ]
        if people[ each_person][ 'mother' ] is None :# if its None so its true 
           
            #probability distribution
            prob_of_gene *=gen_num_prob*trait_prob
        else:  # do this part if false (its not Empty)
            #then use information for parents 
            mami=people[ each_person ][ 'mother']# each percon has one mother 
            father=people[ each_person ][ 'father' ] # each percon has one father
            percent={} 
            for j in [mami, father]: # j to loop over
                num=1 if j in one_gene else 2 if j in two_genes else 0
                perc=0+PROBS[ 'mutation' ] if num==0 else 0.5 if num== 1 else 1 - PROBS['mutation']
                percent [ j ]=perc
            if gen_num==0 :# if 0, that means no one geve a gene
                prob_of_gene *=( 1- percent[ mami ])*( 1- percent[ father ] )
            elif gen_num==1 : # if first if false check this part 
                # 1, that means one of the parents geve gene
                prob_of_gene *=(1-percent[ mami ] )*percent[ father ]+percent[ mami ]*(1- percent[ father] )
            else: # if both contions are false it must be to part
                # 2, that means both  parents geve gene
                prob_of_gene*=percent[ mami ]*percent[ father ] 
            prob_of_gene*=trait_prob
    return prob_of_gene


    raise NotImplementedError


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for pe in probabilities :# pe for person 
        if pe in one_gene: # if the person has one gene then do this part
            gen_num=1# save one into gen
        elif pe in two_genes:# if the person has to ganes then do this part
            gen_num=2# save two in gene
        else :# else here means there no gane 
            gen_num=0# save 0  gene
        # each loop add the vale of probability  
        probabilities[ pe ][ "gene" ][ gen_num ]+=p
        probabilities[ pe ][ "trait" ][ pe in have_trait ]+=p



def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    normalized_rel_prop=probabilities.copy()
    for p in probabilities :# p for person and will go over all probabiltities send foem caller method
        for t in ['gene' , 'trait']:# t as string to check gane and trait
            summed=sum( probabilities[ p][ t].values() ) # add values of prob 
            for category in probabilities[ p ][ t] :
                #calculate the value - store in temp then div by summed 
                temp=probabilities[p ][t ][category ]
                value=temp/summed
                normalized_rel_prop[ p][ t][category ]=value
    return normalized_rel_prop

    raise NotImplementedError


if __name__ == "__main__":
    main()
