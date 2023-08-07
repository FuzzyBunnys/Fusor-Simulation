# Fusor Simulation
This all started out as a way of applying some of the stuff I was learning in my numerical methods and electromagnetics classes at UVic. It's grown a bit subsequently as I found it interesting in its own right. While there's not a ton of simuation work that's been done on Fusors, I'm aware of [some work](https://mattlilley.com/research/electrostatic-fusion/) by Matthew Lilley investigating the dynamics of the fusor. [Liam David](https://fusor.net/board/viewtopic.php?t=14157) has done some top flight work on the Fusor.net forums using particle simulations. This repository serves as a listing of the work I've done so far. Its been completed in essentially two different languages, so that's how I've decided to organize it!
## Python
This python simulation is two dimensional, it uses a finite difference approximation to calculate the electric potential of a metal electrical grid. The current code implements the electric potential of a venetian blind direct energy conversion system. I used the image in my end of degree technical report. A great help in optimizing the speed of this program was the following wonderful [post](https://mattferraro.dev/posts/poissons-equation) by Matt Ferraro. The code is all contained within the SimulationGaussSeidel.py file available on this github repo. 
![2DVoltagePotentialImage](https://raw.githubusercontent.com/FuzzyBunnys/Fusor-Simulation/main/heatmap.png)
## Julia
These are a number of three dimensional simulations of the electric field generated by different grids. It's directly inspired by the work of [Max Lipton](https://e.math.cornell.edu/people/ml2437/). I rewrote his own open source code in Julia as an exercise in learning last Christmas. These simulations are all written in Pluto notebooks, and they operate as follows. A knot shape is parametrized in three dimensional space by some equation and then the electric field is generated via a quadrature integration scheme. A marching cubes algorithm is then used to find a level surface which is plotted to give the final image. There are three notebooks total at the moment but unfortunately I can't figure out how to upload them in a way that allows for their code to work. To get the full effect I would suggest downloading the files from my GitHub and running them locally using Julia. The three notebooks are included in the appropriately name .jl files. 
* Introduction and Unkot 
* The Trefoil
* The Figure Eight
