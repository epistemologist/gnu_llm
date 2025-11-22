# Notes on LLM Mathematical Ability 
We experiment with the mathematical abilities of LLMs in regards to the function gnu(n) = the number of groups of order n up to isomorphism.

## Relevant Links
### On the gnu function
 - [original StackExchange post that inspired this](https://math.stackexchange.com/questions/4931487/does-the-sequence-1-2-3-4-5-6-appear-in-the-number-of-groups-of-order-n/5057635#5057635): does there exist k such that gnu(k+i) = i for i=1..N?
    - relevant OEIS links: [N=2](https://oeis.org/A373648), [N=3](https://oeis.org/A373649), [N=4](https://oeis.org/A373650), [N=5](https://oeis.org/A381335), original post mentions that k = 29436121 gives a solution for i=1..6, no solutions for i=1..8: 
    > I’ve not found a 1,2,3,4,5,6,7 yet - it’s possible 1,2,3,4,5,6 is the longest we can get.
 - conditions on whether gnu(n) = k for small k
    - [Three-group numbers by Olsson](https://web.math.ku.dk/~olsson/manus/three-group-numbers.pdf) gives conditions for k=1,2,3 
    - [Counting groups: gnus, moas and other exotica by Conway et al](https://www.math.auckland.ac.nz/%7Eobrien/research/gnu.pdf) gives conditions for k=1,2,3,4 (see Thm 7.1-4) as well as general discussion of the gnu function 
    - [Orders for Which There Exist Exactly Four or Five Groups by Miller](https://www.pnas.org/doi/abs/10.1073/pnas.18.7.511) gives conditions for k=4,5; written in 1932 and does not use standard notation
    - [Orders for which there exist exactly six or seven groups by Mahmoud](https://arxiv.org/abs/2405.04794) gives conditions for k=6,7
### Calculation of gnu for general n
 - [`NrSmallGroups` in SmallGrp GAP package](https://docs.gap-system.org/pkg/smallgrp/doc/chap1.html#X7C587F2A82BEAD19) for small n, also available in Sage
 - [this GAP package](https://github.com/olexandr-konovalov/gnu) provides algorithms for the general calculation of gnu(n)
### Use of LLMs for Mathematical Conjectures
 - [original AlphaEvolve paper](https://arxiv.org/abs/2511.02864) demonstrating LLMs generating mathematical constructions
 - [OpenEvolve](https://github.com/algorithmicsuperintelligence/openevolve): open source implementation of AlphaEvolve 

## LLM Experiments
 - [k=4, GPT5](https://chatgpt.com/share/691e990c-3c30-8012-a5a9-122ef4b21600): provided with Miller reference, fails to solve
 - [k=1,2,3; gemini-2.5-pro](https://gist.github.com/epistemologist/3aefbcfa41021a16b0e1bcdab0f158b7): when provided with "Three-group numbers" as pdf, solves in one shot
 

