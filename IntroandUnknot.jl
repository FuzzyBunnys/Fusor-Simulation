### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 6eeb2f90-43c4-11eb-2fb0-23fb34f03217
md"## Introduction"

# ╔═╡ 0848c260-43c0-11eb-218b-59a047c6233f
md"This notebook implements [Max Lipton's](https://e.math.cornell.edu/people/ml2437/) Python code on knot electric potentials in [Julia](https://julialang.org/). The original code it is based on can be viewed at the following [link](https://github.com/ml2437/knot-potential-surfaces/blob/master/Electrostatic%20Knot%20Theory.ipynb). Max has also written two papers on this topic, while a lot of the technical details are advanced and I don't understand it all, there's a lot of useful information in them. They deal with the lower bound on the critical points in a [charged knot] (https://e.math.cornell.edu/people/ml2437/tunnel_number_bound.pdf) and the relationship between the [critical sets and Morse code](https://e.math.cornell.edu/people/ml2437/Morse_Code.pdf) of a charged knot. Future plans for this notebook include restating the math in a geometric algebra formalism, implementing posits instead of floats, and optimizing the code for parallel computation. In general I interleave the translated code with markdown comments that explain the math for my own understanding."

# ╔═╡ 36285200-6fe8-11eb-1c4e-dfa1a6692b8a
md"The end goal of all of this work is to have a robust mathematical simulation of a Fusor, a device which I hope will allow me to achieve a self-sustaining fusion reaction."

# ╔═╡ 08b4de70-5165-11eb-01fe-d3732e329757
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.Registry.update()
end

# ╔═╡ 7bd300a0-43cb-11eb-153a-1fe35dc585bb
#add all the packages I want using Pkg.add() then using pkgname
#think I want to use PlutoUI - this gives a bunch of user interface options like 
#sliders and stuff that let you interact with your plots and data.
begin
	Pkg.add("PlutoUI")
	using PlutoUI
end

# ╔═╡ 004c0740-44a4-11eb-2ec2-15a6bd0685a3
begin
	Pkg.add("Plots")
	using Plots
	Pkg.add("PlotlyJS")
	plotly()
end

# ╔═╡ 59570b70-4b91-11eb-252a-db35d8b6e956
begin
	Pkg.add("GeometryBasics")
	Pkg.add("Images")
	using GeometryBasics
	using LinearAlgebra
	using Images
end


# ╔═╡ b5a0612e-4419-11eb-118b-9db02b2086ea
#loading the QuadGK package
begin
	Pkg.add("QuadGK")
	using QuadGK
end

# ╔═╡ 2fc65542-4ae7-11eb-1e6e-5f53acc8d0df
# this might require some additional packages, in particular geometrybasics, geometrytypes and staticarrays
begin
	Pkg.add("Meshing")
	using Meshing
end

# ╔═╡ 50acd190-43ca-11eb-34ad-2bed1a94dd29
md"## Knot Parametrizations"

# ╔═╡ 90400640-43c7-11eb-10a4-831d39da8bde
#Order of accuracy for Gaussian Quadrature - don't know why it shows sometimes 
quadorder = 20 #like python Julia doesn't require ';'

# ╔═╡ 897f19d0-43d7-11eb-392c-6fee25b1eb92
#Discretization of the domain interval for knot parametrizations
t = range(0, stop = 2*pi, length = 1000)# an array of a thousand components

# ╔═╡ 29152b00-43ca-11eb-0221-cd0ccf8aa735
md"We begin by defining the paramaterizations of a number of different knots. In doing so we use the discretization variable t we defined above which is a linear array of 1000 elements from 0 to 2$\pi$. This means that we've swept out an entire circle in 1000 discrete points. The quadorder variable is used for integration later on. We are going to define the following knot types, the Unknot, the Trefoil, the Figure 8, the (3,1) Torus Knot, the (5,1) Torus Knot, the Cinquefoil Knot, the 3-Twist Knot, the Granny Knot, the Square Knot and finally the Endless Knot. Each knot has 3 functions for returning x,y and z coordinates and a numerator function that returns the  derivative. "

# ╔═╡ dd692ae2-440e-11eb-1c97-dd1f4c7e6178
md"#### The Unknot"

# ╔═╡ 061fac90-43db-11eb-3b00-bd5210adbe6e
md"The formula for a circle centered at the origin is of course $r=1$ in polar coordinates. Here to output $x$ and $y$ values we use $x = \cos(t)$ and $y = \sin(t)$. Our value for $z$ is $0$ since the unknot sits on an axis and has no depth."

# ╔═╡ 20492eb0-440e-11eb-06de-c147e3f93660
function unknotnum(t) #the numerator in the integral below
	1
end

# ╔═╡ 5fe42430-440e-11eb-228b-53a81b657786
function unknotx(t)
	cos.(t) #need the dot for elementwise 
end

# ╔═╡ 68cb0aa0-440e-11eb-1565-c94c7e099a0a
function unknoty(t)
	sin.(t)
end

# ╔═╡ 6eba5bf0-440e-11eb-03d6-4d2f8e8fc9f1
function unknotz(t)
	0
end

# ╔═╡ 2b1c2810-449a-11eb-1ea3-2bed99095cff
md"The Unknot is a simple circle, we can visualise the structure using our functions above. "

# ╔═╡ bdd1fb60-493d-11eb-302d-aba1d7e3d74f
plot(unknotx(t),unknoty(t), title = "The Unknot", label = "Unknot", lw = 3, aspect_ratio=:equal)
# want to improve the presentation, ideally it'd be in 3d and scaled symetrically

# ╔═╡ ba86cf00-440e-11eb-033c-8fd03e12333b
md"Having established a parameterization we can now move to integration to compute the electric potential at a specific point (a,b,c). The core task here is how to integrate the electric potential. Recall that the formula for computing an electric potential $\phi_{E}= -\int_{C} E \cdot dl$, which is a line integral. We can ignore the negative for convenience purposes. The electric potential is the voltage at a point in space. As in an electric circuit where we have to measure voltage relative to a reference point (there called ground), here we measure voltage relative to a point at $\infty$ where $\phi=0$. Of course to calculate the electric potential or voltage at a point, we have to know what the electric field is at that point. Having parametrized our knot shapes, we can calculate the electric field by measuring the distance of a point in space from the knot, and then assigning a fixed charge to the knot itself. The Unknot is simply a ring of charge, a special case of which I solved during my Electrical Engineering degree. Now we need to solve the general case in 3 dimensions for the electric field all round the ring of charge rather than at a specific point. "

# ╔═╡ 669c1870-44d9-11eb-09b5-d73398e82adb
md"We will simply state that a point $x$ at some distance $R$ from an electric charge $k$ has an electric potential proportional to $R^{-1}$. This gives us the integral below, where $r(t)$ is the function parametrizing our knot and $x$ is some point. This is the general formula for any knot parametrization, we will need to implement a specific version for the Unknot. "

# ╔═╡ 43aac0b0-4492-11eb-2b05-c51da990da76
md"$\phi=\int_0^{2\pi} \frac{|r'(t)|}{|x-r(t)|}dt$"

# ╔═╡ 6f51d500-4492-11eb-01a6-d9925465c853
md" In our case we're working in $\mathbb{R}^3$ so that means we need to implement a distance including $x,y$ and $z$ coordinates. Our formula for the Unknot will thus end up as below. With the implementation of the function to be integrated written in Julia immediately below it."

# ╔═╡ bc15c6c0-45ab-11eb-2890-ad4730c01e57
md"$\phi=\int_0^{2\pi} \frac{1}{\sqrt{(x-\cos(t))^2+(y-\sin(t))^2+z^2)}}dt$"

# ╔═╡ a33a0e30-4895-11eb-18c8-d3cf673935d9
function unknotintegralfunction(a,b,c,t)
	1/(sqrt(((a.-cos.(t))^2).+((b.-sin.(t))^2).+c^2)) #using .math means element wise
end

# ╔═╡ 17ad58e0-45ac-11eb-31d5-5738c05d21c0
md"We've skipped some details in the derivation, but via cursory examination it should be apparent that we're calculating a distance in $\mathbb{R}^3$. ie: see that we're calculating our distance in the $x$ coordinate by subtracting $cos(t)$, the position of the knot, from our position in $x$. So the denominator of our fraction is the tried and true pythagorean formula. Since the Unknot is centered at the origin and has no height in $z$ we don't subtract anything from it since one's location in $z$ is always the distance from $0$. It's also important to note that the potential we have here is a single number, this derivation doesn't break apart the information into a coordinate based representation."

# ╔═╡ 587dcdd0-440f-11eb-2db6-19abd03ce8c2
md"To actually integrate the formula above we will need an integration function for a single variable in Julia. The best current path seems to be a quadrature routine using the [QuadGK package](https://juliahub.com/ui/Packages/QuadGK/hq5Ol/2.4.1). "

# ╔═╡ 68d75000-4619-11eb-33b2-094cac9629bc
md"There's all sorts of fancy terminology and I need to understand how this works. To restate the basics, I'm integrating a value at a point in 3 dimensions $x$, $y$, and $z$ over a range from $[0,2\pi)$. The single variable is the $t$ of my parametrization function. So in my function below I pass the position (a,b,c) to the function which provides the point in space $t$ over the range specified by the second and third variables. Note the $t \rightarrow$, this tells Julia that the range for quadgk $0$ to $2\pi$ is assigned to the variable $t$. "

# ╔═╡ ded6aae0-4be0-11eb-266c-93a8daeaf540
#setup our storage matrix
A = zeros(41,41,41)

# ╔═╡ 50965470-488d-11eb-3b0a-313cbaf71261
function unknotpotential(a,b,c,A) 
	ans, err = quadgk(t->unknotintegralfunction(a,b,c,t),0,(2*pi),order=quadorder) 
	# we'll use the default errors to start
	A[floor(Int,(10*a+21)),floor(Int,(10*b+21)),floor(Int,(10*c+21))]=ans# need a way to store a,b,c with the ans
	return [a,b,c,ans,err]
end

# ╔═╡ 02920c30-4be2-11eb-244f-4de20ca35dd6
md"We have to generate a 3 dimensional array to store all of the data for our marching cubes algorithm later on. This array will be passed to the functions in question and they'll write the potential value to the specific location in the array. "

# ╔═╡ ecf19c60-4984-11eb-2919-771158df41e2
md"The integration function currently returns the correct value for the potential at the origin. The origin for the unknot is a critical point, so its natural for us to ask why did the integral return the value of $2\pi$ rather than $0$? Recall the definition of a [critical point](https://en.wikipedia.org/wiki/Critical_point_(mathematics)): for a real valued function a critical point is a point where the function is either not differentiable or its derivative is equal to $0$. Above I've calculated the critical value, which is the value of the function at the critical point. So there is a difference between the electric potential, $\phi$ and the electric field, $E = -\nabla\phi$. The electric field is the gradient of the electric potential, so the critical points of the potential are where the **electric field** is $0$, not where the electric potential is $0$. Since the electric potential is a voltage, which by definition is a relative value, the function is saying that there are $2\pi$ volts at the origin relative to a location at $\infty$. "

# ╔═╡ a511cad0-499a-11eb-0c87-e73b17b41f84
md"We can now try and calculate the electric potential for a whole region of space. Unlike older languages such as Matlab or Python you don't need to use a meshgrid type function in Julia. We can instead create our grid in $\mathbb{R}^3$ using list comprehension. This is a fun little programming technique that aims to copy set builder notation. So here I pass a function and a dimension to my inflate function, which then goes and gives all of these values to my passed function. A little bit of Julia syntax to take note of, for the range for $x$,$y$ and $z$ below we have $-1*N:0.1:N$. Which tells us to generate a range from $-1*N$ to $N$ in steps of $0.1$. So if I choose $N=4$ we would have a range $-4:0.1:4$ which would actually be $-4,-3.9,-3.8,...3.8,3.9,4$. We have to be careful with the range we feed this function, since we have a fraction of the form  $\frac{1}{a}$ where $a$ is our distance from the ring, for very small values of $a$ or when $a=0$ we're going to have a problem of blowup. Therefore we need to disallow some values from being included in the set. To do that we can apply some conditions to our list comprehension. In set builder notation we'll write: $\{x,y,z | x,y,z \in \mathbb{R} \wedge x \neq \cos(t) \wedge y \neq \sin(t) \wedge z \neq 0 \}$. Which translates quite nicely and directly into the Julia code one sees below. "

# ╔═╡ 4a56c450-4a0e-11eb-2432-fde548b9f590
#we're going to plagiarize this function from Tomas' answer in google groups
#N should be a positive integer
inflateunknot(f,N,p,A)=[f(x,y,z,A) for x in -1*N:p:N, y in -1*N:p:N, z in -1*N:p:N if z ≠ 0 || (x*x +y*y)≠ 1 ]

# ╔═╡ e7a1b8be-4a1b-11eb-16cc-0b68c71c2418
#now we're going to call the function itself
ans = inflateunknot(unknotpotential,2,0.1,A)

# ╔═╡ 7030e220-4a6b-11eb-1db0-6399f284f4bc
md"Having successfully calculated the electric potential at all of these points we now have an enormous tuple. The tuple is itself made up of tuples. The first three elements of each element in the tuple are the position in space and the last two are the electric potential and the error in the integration calculation.  "

# ╔═╡ 379a5fb0-4be5-11eb-1338-f360ec868a99
A #Just confirming data was written to the array

# ╔═╡ 5e06cf02-4a6c-11eb-24dd-6b1bceb4f832
md"The next step is to visualize the data. Max calculated a level potential surface. Recall the definition of a [level surface](https://web.ma.utexas.edu/users/m408m/Display12-6-4.shtml): For a function $\phi = f(x,y,z): U \subseteq \mathbb{R}^3 \rightarrow \mathbb{R}$ the level surface of value $c$ is the surface $s$ in $U \subseteq \mathbb{R}^3$ on which $\phi \rvert_{s}=c$. So in our case we're collecting all the points that have the same voltage."

# ╔═╡ 89db1b90-4a71-11eb-0280-e9986755c8ef
md"The trick of course is to generate the surface to plot. How do we collect those points? To do that Max used the [marching cubes algorithm](https://en.wikipedia.org/wiki/Marching_cubes). There is of course a package that gives us this capability in Julia called [meshing](https://juliageometry.github.io/Meshing.jl/stable/). There's a lovely introduction to marching cubes and the isosurface problem at the following [link](https://0fps.net/2012/07/12/smooth-voxel-terrain-part-2/), and if you'd like a very good exposition of the technical details involved you can investigate the classic explanation by [Paul Bourke](http://paulbourke.net/geometry/polygonise/). Very basically the marching cubes algorithm takes some sort of input, either an equation or an array, and outputs a mesh. The technical name for this problem is isosurface extraction."

# ╔═╡ b19d20b0-4b25-11eb-2d21-e55d46ad40ee
md"The first bit of important information for us is to note that the scalar values that we're going to use to generate our manifold are stored at each vertex. Recall from graph theory that a vertex is a node, and an edge is just an edge. So for a generic cube, we'd have 8 vertices (corners), and 12 edges connecting all the vertices. What the marching cubes algorithm does is generate a surface based on a boundary the user specifies. So I should be able to specify a value for the electric potential, and find a boundary around the ring where that value holds. To do that we need to grab the matrix A we defined earlier and feed it into the marching cubes algorithm. "

# ╔═╡ a4cac4a0-4be8-11eb-04e6-a7328783b6aa
md"We'll just copy what Max did for the level value of our isosurface, which is to take the potential value at the origin and arbitrarily add 0.5 to it." 

# ╔═╡ 5b0caa10-4be6-11eb-2047-6deab2f3b7fd
c=A[20,20,20]+0.5

# ╔═╡ ba45dac0-4b27-11eb-1675-ef9d2ae9f9ca
#generate a mesh
points,surfaces = isosurface(A,MarchingCubes(iso=c), origin = Vec(1,1,1), widths = Vec(4.0,4.0,4.0)) 

# ╔═╡ d4f6ebe0-4be8-11eb-1239-3320defd23ba
md"Now that the points and surfaces of the mesh have been generated we need to plot them! To do that we can use the Mesh3d function. And to access the values in my points tuple I'm going to have to use a list comprehension. Having provided the points mesh3d will autogenerate the mesh for me. There are currently some issues getting mesh3d to run with the plotting backend I'm using. PlotlyJS is not currently setup for Pluto so until that happens I think this will have some issues. In addition for the kind of shape I'm generating I have to change the alphahull parameter in mesh3d. It should be set to 5." 

# ╔═╡ f215b4d0-4bfd-11eb-0a33-d5711fbef747
mesh3d([a for (a,_,_) in points], [b for (_,b,_) in points], [c for (_,_,c) in points],alphahull=5)

# ╔═╡ a8d9d032-4d88-11eb-06f8-01374ce6b505
md"A quick and dirty way to display the data is with a scatter plot. Again a list comprehension is used to access the individual elements in the tuple points."

# ╔═╡ f4917832-4d87-11eb-0f91-532328b1a8f6
scatter([a for (a,_,_) in points],[b for (_,b,_) in points], [c for (_,_,c) in points] , title = "Unknot Potential", ms = 1)

# ╔═╡ b2bc6ea0-4d2e-11eb-32da-bb24b71c18c4
md"This show method along with the BWImage function were copied from the following [spot](https://gist.github.com/pbouffard/3d48d3c47d9bd70e7c9f52f984d14245). It apparently comes from the [2020 JuliaCon talk](https://www.youtube.com/watch?v=IAF8DjrQSSk). Either way now that it's included I've got a very quick and dirty way to view a 2D matrix as an image." 

# ╔═╡ 522123b0-4d2e-11eb-3fc4-c946ba08cbaa
begin
	struct BWImage
		data::Array{UInt8, 2}
		zoom::Int
	end
	function BWImage(data::Array{T, 2}; zoom::Int=1) where T <: Real
		BWImage(floor.(UInt8, clamp.(((data .- minimum(data)) / (maximum(data) .- minimum(data))) * 255, 0, 255)), zoom)
	end
	
	import Base: show
	
	function show(io::IO, ::MIME"image/bmp", i::BWImage)

		orig_height, orig_width = size(i.data)
		height, width = (orig_height, orig_width) .* i.zoom
		datawidth = Integer(ceil(width / 4)) * 4

		bmp_header_size = 14
		dib_header_size = 40
		palette_size = 256 * 4
		data_size = datawidth * height * 1

		# BMP header
		write(io, 0x42, 0x4d)
		write(io, UInt32(bmp_header_size + dib_header_size + palette_size + data_size))
		write(io, 0x00, 0x00)
		write(io, 0x00, 0x00)
		write(io, UInt32(bmp_header_size + dib_header_size + palette_size))

		# DIB header
		write(io, UInt32(dib_header_size))
		write(io, Int32(width))
		write(io, Int32(-height))
		write(io, UInt16(1))
		write(io, UInt16(8))
		write(io, UInt32(0))
		write(io, UInt32(0))
		write(io, 0x12, 0x0b, 0x00, 0x00)
		write(io, 0x12, 0x0b, 0x00, 0x00)
		write(io, UInt32(0))
		write(io, UInt32(0))

		# color palette
		write(io, [[x, x, x, 0x00] for x in UInt8.(0:255)]...)

		# data
		padding = fill(0x00, datawidth - width)
		for y in 1:orig_height
			for z in 1:i.zoom
				line = vcat(fill.(i.data[y,:], (i.zoom,))...)
				write(io, line, padding)
			end
		end
	end
end

# ╔═╡ 5bed8820-4d2e-11eb-0e49-51285dad145e
# A black and white image of the potential
BWImage(A[:,20,:], zoom=10)

# ╔═╡ 6b5d40c0-4d2e-11eb-3ffc-3b2857b90f7b
#A different view of the potential
BWImage(A[:,:,20], zoom=10)

# ╔═╡ 2fa8fa52-4d39-11eb-1014-4f20761386d3
BWImage(A[20,:,:], zoom=10)

# ╔═╡ e158b030-48d2-11eb-3a88-5f96538a8c9c
md"Following the calculation of the actual electric potential, I want to analyze the orbits of the vector field. The classic equation is of course $F=ma$ for the forces on a particle changing its acceleration. I've calculated the electric potential above. In this case however I obviously have to take into account the charge on the particle itself. Since I'm not considering magnetic effects at the moment I have a nice simple equation $F=qE$, where $E=-\nabla \phi$. Putting that all together we'd have $ma = -q\nabla\phi$. So first I need to get a 3d matrix that represents the electric field, at the same time I want to express my equations symbolically so its a bit easier to reason with them. To do that I need to setup a number of different equations for the $x,y$ and $z$ components of the force the charged particle will feel. Here's our earlier equations repurposed for finding $E$."

# ╔═╡ 8babc330-6fdc-11eb-2591-0f503e7efae4
md"$E=-\nabla \phi=-\nabla \int_0^{2\pi} \frac{|r'(t)|}{|x-r(t)|}dt=-\nabla \int_0^{2\pi} \frac{1}{\sqrt{(x-\cos(t))^2+(y-\sin(t))^2+z^2)}}dt$"

# ╔═╡ bb5900fe-6fdd-11eb-1970-eb8ce136fa42
md"Since we're dealing with the original use case for vector calculus, and this has all been checked and confirmed by lots of mathematicians we can go ahead and the integral with the gradient operator. My use of notation here is imprecise for sure but I think it suffices to get the general point across. So I'm going to have three equations where I'm calculating $E$ in $x,y$ and $z$. Since the electric field depends on the distance from the charge source, essentially what we have here is a series of distance equations in a single direction of interest. So the $x$ component of the electric field is how far our position in $x$ is from our circle of charge. "

# ╔═╡ 0345b880-4934-11eb-326d-79484684c8e9
md"$E_x = -\frac{d\phi}{dx}=-\frac{1}{x-\cos(t)}$"

# ╔═╡ 938949f0-701e-11eb-11b1-f77c4c98bdd1
function unknotefieldx(x,t)
end

# ╔═╡ 11ad7700-4934-11eb-1e5d-e79001f7c099
md"$E_y=-\frac{d\phi}{dy}=-\frac{1}{y-\sin(t)}$"

# ╔═╡ 94a9ae0e-701e-11eb-1dab-cf272740221f
function unknotefieldy(y,t)
end

# ╔═╡ 12abad20-4934-11eb-1242-19509302b73e
md"$E_z=-\frac{d\phi}{dz}=-\frac{1}{z}$"

# ╔═╡ 95689140-701e-11eb-0cdf-0b0ebf9f48bc
function unknotefieldz(z)
	-1/z
end

# ╔═╡ 4459aab0-493a-11eb-1982-1965fec7b47e
md"Now that we've got these 3 equations setup we can start searching for some zeros! Since Max has already done the work we know that the only critical point for this potential is at the origin. Is it possible to show that algebraically using the techniques I learned in my ODE's course?"

# ╔═╡ f09e3bf0-511f-11eb-3484-59558508d5a2
md"By observation we can see that the equation for $d\phi/dz$ is independent of the the other two and $d\phi/dx$ and $d\phi/dy$ share a common parameter in $w$. If we set $z=0$ we can see we have a point where our function is undefined. I believe this meets the definition of not differentiable, so we know our critical point has to lie in that plane. Max's [paper](https://e.math.cornell.edu/people/ml2437/tunnel_number_bound.pdf) sets a lower bound on the number of critical points of an electrical field of a knot. I'm not going to break down his proof, since I don't understand it at the moment, but I am going to note some facts about it. The [tunneling number](https://en.wikipedia.org/wiki/Tunnel_number) of a knot, which is an invariant 'fact' about a knot. In the case of the unknot it's tunneling number is $0$. If we take Max's formula $2t(K)+2$ for the minimum number of critical points we find $t=0$ for the unknot and obtain $2$ as the number of critical points for the unknot. We should note that this number of critical points includes the point at $\infty$ since by definition we set $\infty=0$. Thus we have one point at $\infty$ and another somewhere in the plane where $z=0$. Note that $t(K)$ is a function where $K$ is the knot type and $t$ is the tunneling number associate with that knot. "

# ╔═╡ 24270440-5137-11eb-1882-79b5158aa440
md"So how do we determine where in the plane $z=0$ the critical point is itself located? I know that $t$ is our variable of integration from $0 \rightarrow 2\pi$. And I know from examination of a plot of $sin(t)$ and $cos(t)$ that the total area of these two functions is $0$ for a complete period from $0 \rightarrow 2\pi$. So if I choose $x=0$ and $y=0$ well then I have $\frac{1}{0-\sin(t)}$ and $\frac{1}{0-\cos(t)}$ for all $t$. I just put these together and I should see intuitively why the origin is the critical point."

# ╔═╡ 643962f0-5153-11eb-1ca8-7935796dce93
md"Now what is the class of the critical point at the center of the unknot? I believe I need to calculate the Jacobian matrix now and examine it's eigenvalues. Recall the definition of a Jacobian: "

# ╔═╡ 14393ad0-515f-11eb-32f8-1949428e445d
md"$\mathbf{J}=\begin{bmatrix}
E_{xx}&E_{xy}&E_{xz}\\
E_{yx}&E_{yy}&E_{yz}\\
E_{zx}&E_{zy}&E_{zz}\\
\end{bmatrix}$"

# ╔═╡ 31c386e0-5160-11eb-095e-a59089a4a0df
md"Which for us gives us the following matrix. I think this is correct, I'd need to confirm with someone who knows a bit more about this than me though..."

# ╔═╡ d3221020-515f-11eb-3927-81bb6fa164c2
md"$\mathbf{J}=\begin{bmatrix}
-\frac{1}{(x-\cos(t))^2}&0&0\\
0&-\frac{1}{(y-\sin(t))^2}&0\\
0&0&-\frac{1}{z^2}\\
\end{bmatrix}$"

# ╔═╡ 3d9b72c0-5160-11eb-0e16-03e529082f25
md"Since this is a diagonal matrix we can simply read the eigenvalues off the diagonal. Now we have to evaluate them at our critical point $(0,0,0)$."

# ╔═╡ 864c6310-5162-11eb-15ea-613e8922a9c4
md"$\lambda_1=-\frac{1}{(x-\cos(t))^2}=-\frac{1}{\cos^2(t)}$"

# ╔═╡ 86b9f1f0-5162-11eb-2ec3-9d73a6168aca
md"$\lambda_2=-\frac{1}{(y-\sin(t))^2}=-\frac{1}{\sin^2(t)}$"

# ╔═╡ 8761f1c2-5162-11eb-3320-4f1a9e5b829f
md"$\lambda_3=-\frac{1}{z^2}=-\frac{1}{0}$"

# ╔═╡ f2792590-5163-11eb-3308-c791a252c7be
md"So now how do I evaluate these results? I know I can use the techniques, for negative real, complex etc. But what category to my points fall into? I have three results, and if I evaluate for the range of $t$ for my knot then the values are $\lambda_1 = -\frac{1}{\pi}$,$\lambda_2 = -\frac{1}{\pi}$ and $\lambda_3 = -\frac{1}{0}$. So I have two real equal negative eigenvalues and a third negative real eigenvalue (I think, not 100% sure about that last fraction). According to Rod, if all my values are negative and real then the point is asymptotically stable in $R^n$. "

# ╔═╡ b18cd05e-6fed-11eb-1962-253005ff4079
md"So we can now generate the actual electric field to see what happens with our particle. Since I have the equations for the electric field based on it's position in space, I can just do a leapfrog scheme to figure out how it moves. I don't have to have a grid that stores my electric field. I can feed an initial position to the particle, increment that position and then check it's new position. "

# ╔═╡ cacae870-4def-11eb-195a-9944787e7f43
md"Rather than continue extending this already rather long notebook, I'm going to generate a different notebook for each type of knot and work with them there."

# ╔═╡ Cell order:
# ╟─6eeb2f90-43c4-11eb-2fb0-23fb34f03217
# ╟─0848c260-43c0-11eb-218b-59a047c6233f
# ╟─36285200-6fe8-11eb-1c4e-dfa1a6692b8a
# ╠═08b4de70-5165-11eb-01fe-d3732e329757
# ╠═7bd300a0-43cb-11eb-153a-1fe35dc585bb
# ╠═004c0740-44a4-11eb-2ec2-15a6bd0685a3
# ╠═59570b70-4b91-11eb-252a-db35d8b6e956
# ╟─50acd190-43ca-11eb-34ad-2bed1a94dd29
# ╠═90400640-43c7-11eb-10a4-831d39da8bde
# ╠═897f19d0-43d7-11eb-392c-6fee25b1eb92
# ╟─29152b00-43ca-11eb-0221-cd0ccf8aa735
# ╟─dd692ae2-440e-11eb-1c97-dd1f4c7e6178
# ╠═061fac90-43db-11eb-3b00-bd5210adbe6e
# ╠═20492eb0-440e-11eb-06de-c147e3f93660
# ╠═5fe42430-440e-11eb-228b-53a81b657786
# ╠═68cb0aa0-440e-11eb-1565-c94c7e099a0a
# ╠═6eba5bf0-440e-11eb-03d6-4d2f8e8fc9f1
# ╟─2b1c2810-449a-11eb-1ea3-2bed99095cff
# ╠═bdd1fb60-493d-11eb-302d-aba1d7e3d74f
# ╟─ba86cf00-440e-11eb-033c-8fd03e12333b
# ╟─669c1870-44d9-11eb-09b5-d73398e82adb
# ╟─43aac0b0-4492-11eb-2b05-c51da990da76
# ╟─6f51d500-4492-11eb-01a6-d9925465c853
# ╟─bc15c6c0-45ab-11eb-2890-ad4730c01e57
# ╠═a33a0e30-4895-11eb-18c8-d3cf673935d9
# ╟─17ad58e0-45ac-11eb-31d5-5738c05d21c0
# ╟─587dcdd0-440f-11eb-2db6-19abd03ce8c2
# ╠═b5a0612e-4419-11eb-118b-9db02b2086ea
# ╟─68d75000-4619-11eb-33b2-094cac9629bc
# ╠═ded6aae0-4be0-11eb-266c-93a8daeaf540
# ╠═50965470-488d-11eb-3b0a-313cbaf71261
# ╟─02920c30-4be2-11eb-244f-4de20ca35dd6
# ╟─ecf19c60-4984-11eb-2919-771158df41e2
# ╟─a511cad0-499a-11eb-0c87-e73b17b41f84
# ╠═4a56c450-4a0e-11eb-2432-fde548b9f590
# ╠═e7a1b8be-4a1b-11eb-16cc-0b68c71c2418
# ╟─7030e220-4a6b-11eb-1db0-6399f284f4bc
# ╠═379a5fb0-4be5-11eb-1338-f360ec868a99
# ╟─5e06cf02-4a6c-11eb-24dd-6b1bceb4f832
# ╟─89db1b90-4a71-11eb-0280-e9986755c8ef
# ╠═2fc65542-4ae7-11eb-1e6e-5f53acc8d0df
# ╟─b19d20b0-4b25-11eb-2d21-e55d46ad40ee
# ╟─a4cac4a0-4be8-11eb-04e6-a7328783b6aa
# ╠═5b0caa10-4be6-11eb-2047-6deab2f3b7fd
# ╠═ba45dac0-4b27-11eb-1675-ef9d2ae9f9ca
# ╟─d4f6ebe0-4be8-11eb-1239-3320defd23ba
# ╠═f215b4d0-4bfd-11eb-0a33-d5711fbef747
# ╟─a8d9d032-4d88-11eb-06f8-01374ce6b505
# ╠═f4917832-4d87-11eb-0f91-532328b1a8f6
# ╟─b2bc6ea0-4d2e-11eb-32da-bb24b71c18c4
# ╟─522123b0-4d2e-11eb-3fc4-c946ba08cbaa
# ╠═5bed8820-4d2e-11eb-0e49-51285dad145e
# ╠═6b5d40c0-4d2e-11eb-3ffc-3b2857b90f7b
# ╠═2fa8fa52-4d39-11eb-1014-4f20761386d3
# ╟─e158b030-48d2-11eb-3a88-5f96538a8c9c
# ╟─8babc330-6fdc-11eb-2591-0f503e7efae4
# ╟─bb5900fe-6fdd-11eb-1970-eb8ce136fa42
# ╟─0345b880-4934-11eb-326d-79484684c8e9
# ╠═938949f0-701e-11eb-11b1-f77c4c98bdd1
# ╟─11ad7700-4934-11eb-1e5d-e79001f7c099
# ╠═94a9ae0e-701e-11eb-1dab-cf272740221f
# ╟─12abad20-4934-11eb-1242-19509302b73e
# ╠═95689140-701e-11eb-0cdf-0b0ebf9f48bc
# ╟─4459aab0-493a-11eb-1982-1965fec7b47e
# ╟─f09e3bf0-511f-11eb-3484-59558508d5a2
# ╟─24270440-5137-11eb-1882-79b5158aa440
# ╟─643962f0-5153-11eb-1ca8-7935796dce93
# ╟─14393ad0-515f-11eb-32f8-1949428e445d
# ╟─31c386e0-5160-11eb-095e-a59089a4a0df
# ╟─d3221020-515f-11eb-3927-81bb6fa164c2
# ╟─3d9b72c0-5160-11eb-0e16-03e529082f25
# ╟─864c6310-5162-11eb-15ea-613e8922a9c4
# ╟─86b9f1f0-5162-11eb-2ec3-9d73a6168aca
# ╟─8761f1c2-5162-11eb-3320-4f1a9e5b829f
# ╟─f2792590-5163-11eb-3308-c791a252c7be
# ╠═b18cd05e-6fed-11eb-1962-253005ff4079
# ╟─cacae870-4def-11eb-195a-9944787e7f43
