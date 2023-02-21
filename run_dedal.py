"""================================================================================================
Test script for DEDAL model. Place dedal folder in the same directory as this script.

Ben Iovino  02/21/23   VecAligns
================================================================================================"""


import tensorflow as tf
import tensorflow_hub as hub
from dedal import infer


def main():
    """=============================================================================================
    Run the DEDAL model to get a pairwise alignment between two proteins.
    ============================================================================================="""

    # Load model and preprocess inputs
    dedal_model = hub.load('https://tfhub.dev/google/dedal/3')
    protein_a = 'SVCCRDYVRYRLPLRVVKHFYWTS'
    protein_b = 'VKCKCSRKGPKI'
    inputs = infer.preprocess(protein_a, protein_b)
    with open('inputs.txt', 'w', encoding='utf8') as f:
        f.write(str(inputs))

    align_out = dedal_model(inputs)
    with open('align_out.txt', 'w', encoding='utf8') as f:
        f.write(str(align_out))

    output = infer.expand(
        [align_out['sw_scores'], align_out['paths'], align_out['sw_params']])
    output = infer.postprocess(output, len(protein_a), len(protein_b))
    alignment = infer.Alignment(protein_a, protein_b, *output)
    with open('alignment.txt', 'w', encoding='utf8') as f:
        f.write(str(alignment))


if __name__ == '__main__':
    main()
