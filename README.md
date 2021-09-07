# Fusor Simulation
This all started out as a way of applying some of the stuff I was learning in my numerical methods and electromagnetics classes at UVic. It's grown a bit subsequently as I found it interesting in its own right. There's been [some work](https://mattlilley.com/research/electrostatic-fusion/) by other people in an attempt to simulate the dynamics of the fusor, but nobody so far as I'm aware, has a really good mathematical model of what's going on in the device.  The work done here has been completed in essentially two different languages, so that's how I've decided to organize it!
## Python
This python simulation is two dimensional, it uses a finite difference approximation to calculate the electric potential of a metal electrical grid. The current code implements the electric potential of a venetian blind direct energy conversion system. I used the image in my end of degree technical report. A great help in optimizing the speed of this program was the following wonderful [post](https://mattferraro.dev/posts/poissons-equation) by Matt Ferraro. The code is all contained within the SimulationGaussSeidel.py file available on this github repo. 
![2DVoltagePotentialImage](https://raw.githubusercontent.com/FuzzyBunnys/Fusor-Simulation/main/heatmap.png)
## Julia
These are a number of three dimensional simulations of the electric field generated by different grids. It's directly inspired by the work of [Max Lipton](https://e.math.cornell.edu/people/ml2437/). These simulations are all written in Pluto notebooks, and if I've configured things correctly should be visible to you, the viewer. The explanations on the notebooks are more detailed, but they operate as follows. A knot shape is parametrized in three dimensional space by some equation and then the electric field is generated via a quadrature integration scheme. A marching cubes algorithm is then used to find a level surface which is plotted to give the final image. 
