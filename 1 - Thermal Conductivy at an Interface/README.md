# Problem 1

## Thermal Conductivity Evaluation at an Interface

Please open the [PDF](https://github.com/ggwadera/cfd-matlab-problems/blob/master/1%20-%20Thermal%20Conductivy%20at%20an%20Interface/Problem%201.pdf) included in this folder to understand about the implementation and results.

The proposed problem consists of a wall composed by two materials
of different thermal conductivities, with a prescribed heat flow *q*'' = 6000 W/m<sup>2</sup>
on the left border, and a convective heat flow on the right border, with *h* = 100 W/(m<sup>2</sup> K) and *T*<sub>∞</sub> = 40 °C. For this case, steady state is assumed.

![Problem representation.](figures/problema.svg)

It is requested to find the temperature profile along the wall using two different methods to evalute the thermal conductivity at the materials interface: thermal resistances, and linear variation. It's also requested to compare the results of both methods against the exact analytical solution.
