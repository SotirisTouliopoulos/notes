# Notes on the steady-state assumption in metabolic networks


## The steady-states assumption

A steady state is an equilibrial condition in a dynamic process in which the rate of input of a state variable is equal to the rate of output of that variable.
For a metabolic pathway at steady state, the rate of input into the pathway, the rate of conversion of A to B and the rate of output are all equal

When we analyze a metabolic network and obtain reactions fluxes, we get the fluxes in the steady-state condition.

If we change a reaction's flux, then we also alter the concentration of other metabolites in the system.
These metabolites must find a way to return to the original steady state or settle to a new one.

Thus, specific pathways are activated that increase or decrease the concentration of the metabolites.


## Correlated reactions

Flux sampling is an unbiased method to analyze the possible states of a metabolic network.
Sampling returns a numeric matrix as output with the flux distribution of each reaction.
A histogram of a flux distribution can be seen in the figure below:

![PFK-optimal-distribution](img/)

Statistical analysis of fluxes distributions can give insights on whether 2 pairwise reactions are correlated or not.
This can be done with a pearson coefficient or with a copula indicator.

Certain reactions have distributions with very low flux values. Their mean value is close to 0 and they have a low variance value too.
An example of this, is the distribution of the `FBP` reaction when sampling under the optimal conditions:

![FBP-optimal-distribution](img/)

Correlation tests to compare this distributions with others are nonsense and they may produce false information.
So we can identify this type of distributions and ingore them when applying correlation tests.

Another method to identify correlation between pairwise reactions is to compute their copula and graphically represent it.
A visual representation of a copula can be seen in the figure below:

![AAA-BBB-copula](img/)

Copulas show the values a reaction can take at different states of its pair reaction.
In the figure above we see that ...

From what is stated, copulas are expected to represent positive, negative or no (when uniform) correlation.
However, in some cases, copulas indicate correlation only when a reaction takes its lowest or highest flux values.
An example of this type of copula is:

![PFK-FBP-copula](img/)

We observe that `PFK` takes its highest values only when `FBP` takes its highest values as well.
This can be explained by the steady-state assumption. `PFK` is a reaction that is involved in the glycolytic pathway,
whereas `FBP` is involved in glyconeogenesis. When `PFK` is highly active it consumes `D-Fructose-6-phosphate` to produce
`D-Fructose-1-6-bisphosphate`. The excessive amount of `D-Fructose-1-6-bisphosphate` satisfies the biomass needs of the cell
and thus, the remaining amount is converted back to `D-Fructose-6-phosphate` through the `FBP` to reach its steady-state concentration.

We can better understand this concept, by exploring the model behavior under sub-optimal conditions. When the cell does not prioritize biomass
production pathways that are opposed to the biomass reaction have increased activity. As such, in this case `FBP` has increased activity,
as the amount of `D-Fructose-1-6-bisphosphate` is increased too. When this happens, `PFK` and `FBP` have increased correlation too.

