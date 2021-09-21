### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 6c1821f0-5234-11eb-369b-d1a90e15cd01
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.Registry.update()
end

# ╔═╡ 0db90400-5237-11eb-1a17-6d996e294d68
begin
	Pkg.add("Plots")
	using Plots
	Pkg.add("PlotlyJS")
	plotly()
end

# ╔═╡ 12ed7810-52ab-11eb-3805-0166f921531c
#loading the QuadGK package
begin
	Pkg.add("QuadGK")
	using QuadGK
end

# ╔═╡ 58612f00-52be-11eb-1376-fd015c444153
begin
	Pkg.add("GeometryBasics")
	Pkg.add("Images")
	using GeometryBasics
	using LinearAlgebra
	using Images
end

# ╔═╡ 6d874ea0-52be-11eb-3357-af01439bf137
begin
	Pkg.add("Meshing")
	using Meshing
end

# ╔═╡ de1e7290-522f-11eb-3c91-8546a8d814bf
md"### The Figure 8"

# ╔═╡ b6b7d1a0-5230-11eb-3fc3-751e3a2dfcfa
md"We start as before with our parametrization of the knot in $x,y$ and $z$. Along with the numerator of our integral $r'(t)$."

# ╔═╡ 772bfa80-5234-11eb-02a3-f185e29c63a9
md"$x=(2+\cos(2t))\cos(3t)$"

# ╔═╡ 7afa67a0-5234-11eb-2fd3-0d8d977b634b
function figure8x(t)
	return (2 .+ cos.(2 .*t)) .*cos.(3 .*t)
end

# ╔═╡ 79d325b0-5234-11eb-3e6c-1324a5470226
md"$y=(2+ \cos(2t))\sin(3t)$"

# ╔═╡ 7b878c70-5234-11eb-39dd-9dc1bbfd76f3
function figure8y(t)
	return (2 .+ cos.(2 .*t)) .*sin.(3 .*t)
end

# ╔═╡ 7a5722c0-5234-11eb-21c9-ef24bfde97b2
md"$z=\sin(4t)$"

# ╔═╡ 7c1cc790-5234-11eb-0c3e-0381bf9a6b02
function figure8z(t)
	return sin.(4 .*t)
end

# ╔═╡ 989fab30-5234-11eb-1360-176b9497415b
md"$r'(t)=40+36\cos(2t)+5\cos^2(2t)+16\cos^2(4t)$"

# ╔═╡ 9758d350-5234-11eb-0470-7de7fb74ebbd
function figure8num(t)
	return 40 .+ 36 .*cos.(2 .*t) .+ 5 .*(cos.(2 .*t))^2 .+ 16 .* (cos.(4 .*t))^2
end

# ╔═╡ f1371bb0-523f-11eb-37cb-c1b9ab9ee329
md"Then we need to define our range for our parameter $t$." 

# ╔═╡ 218e2380-5240-11eb-27e7-374a37c7100f
t = 0:(2*pi/1000):2*pi

# ╔═╡ 43c4d610-5240-11eb-3bbf-993a98d50ea5
md"We then plot our knot."

# ╔═╡ 4dcac110-5240-11eb-2ce5-df407e48f332
plot(figure8x(t),figure8y(t),figure8z(t),title = "Figure 8 Knot")

# ╔═╡ b7c48182-5242-11eb-2eea-7f7090854eaa
md"As usual our order of integration is $5$."

# ╔═╡ b31180c0-5242-11eb-21a0-1175781341f2
quadorder=5

# ╔═╡ 7ac2b430-5294-11eb-180b-c7e8a015731a
md"We may now proceed with defining our integral for the figure 8 knot."

# ╔═╡ 9aa625f0-52aa-11eb-0ed6-eb5a4723e806
md"$\phi=\int_0^{2\pi} \frac{\sqrt{40+36\cos(2t)+5\cos^2(2t)+16\cos^2(4t)}}{\sqrt{(x-(2+\cos(2t))\cos(3t))^2+(y-(2+ \cos(2t))\sin(3t))^2+(z-\sin(4t))^2}}$"

# ╔═╡ 97b5a370-52af-11eb-2029-ff804f96a261
md"Which as before gives us this nasty function with all sorts of trigonometric components, just remember that the top is $|r'(t)|$ and the bottom is $\sqrt{a^2+b^2+c^2}$."

# ╔═╡ 103a8d10-52ab-11eb-14a5-8915384a290b
function figure8integral(a,b,c,t)
return	sqrt(figure8num(t))/sqrt((a-figure8x(t))^2 + (b-figure8y(t))^2 + (c-figure8z(t))^2)
end

# ╔═╡ 5d1c19a0-52ab-11eb-044a-2f0d527709f5
function figure8potential(a,b,c,N,p,A)
	ans, err = quadgk(t->figure8integral(a,b,c,t),0,(2*pi),order=quadorder) 
	# we'll use the default errors to start
	A[floor(Int,((1/p)*a+((N/p)+1))),floor(Int,((1/p)*b+((N/p)+1))),floor(Int,((1/p)*c+((N/p)+1)))]=ans
	return A
end

# ╔═╡ 809e81e0-52b2-11eb-22f5-e12d7e2a10b2
md"We can use the same trick we employed for the trefoil knot to check if a point is on the knot when we integrate. We find the inverse of $z$ and test the outputs against our values for $x$ and $y$."

# ╔═╡ 5c0d2430-52b7-11eb-1a30-ed8986c0a8fb
md"$z^{-1}  = \frac{\sin^{-1}(z)}{4}= t$"

# ╔═╡ 90734ba0-52b7-11eb-0723-652371f58cc1
function testfig8z(x,y,z,t) #so this test function will return true or false
	if z>1
		return true
	elseif z<-1
		return true
	else
		t= (asin(z))/4
		a = figure8x(t)
		b = figure8y(t)
		if x==a && y==b
			return false # the point is on the knot so GET OUTTA HERE
		else
			return true
		end
	end
end

# ╔═╡ a6e01c10-52b7-11eb-1366-05747462a345
function inflatefigure8(f,N,p,t)
	#dynamically setup our storage array
	A = zeros(floor(Int,(2*(N/p)+1)),floor(Int,(2*(N/p)+1)),floor(Int,(2*(N/p)+1)))
	[f(x,y,z,N,p,A) for x in -1*N:p:N, y in -1*N:p:N, z in -1*N:p:N if 	testfig8z(x,y,z,t)]
	return A
end

# ╔═╡ ab5de9c2-52bc-11eb-1186-292123827b72
ans = inflatefigure8(figure8potential,4,0.1,(0:0.1:2*pi))

# ╔═╡ abf1772e-52bc-11eb-33b4-c39608653da6
c = ans[40,40,40]+0.5

# ╔═╡ 03a3f4fe-52bf-11eb-3fe5-1d21d702651d
md"Again we use marching cubes to extract the surface, not 100% sure about the setting of the origin and widths but for now it's working so there's no real need to change it."

# ╔═╡ 55702992-52be-11eb-07b9-b57e098f7b3c
#marching cubes to extract the surface
points,surfaces = isosurface(ans,MarchingCubes(iso=c), ) 

# ╔═╡ 733df1d0-52c0-11eb-2217-43b21eed1187
md"Currently the potentials are doing well, and my initial matrix size  of 41x41x41 was large enough for my unknot, but for both the trefoil and the figure 8 I'm missing a lot of detail in the potentials since I don't have a wide enough field of view. I need to increase the size of A."

# ╔═╡ d6f36a82-52bf-11eb-294b-e5276cabc27b
scatter([a for (a,_,_) in points],[b for (_,b,_) in points], [c for (_,_,c) in points] , title = "Figure 8 Potential", ms = 1)

# ╔═╡ 2bf99490-52e4-11eb-0ddb-2b9c799865e7
md"At the moment everything is working, but from what I can tell I'm generating three copies of my array for storing my data. No clue why that is, it's definitely something to investigate in the future."

# ╔═╡ 4b66edf0-52e4-11eb-25ab-0978549525c6
md"### Critical Points"

# ╔═╡ Cell order:
# ╠═6c1821f0-5234-11eb-369b-d1a90e15cd01
# ╠═0db90400-5237-11eb-1a17-6d996e294d68
# ╠═12ed7810-52ab-11eb-3805-0166f921531c
# ╠═58612f00-52be-11eb-1376-fd015c444153
# ╠═6d874ea0-52be-11eb-3357-af01439bf137
# ╟─de1e7290-522f-11eb-3c91-8546a8d814bf
# ╟─b6b7d1a0-5230-11eb-3fc3-751e3a2dfcfa
# ╟─772bfa80-5234-11eb-02a3-f185e29c63a9
# ╠═7afa67a0-5234-11eb-2fd3-0d8d977b634b
# ╟─79d325b0-5234-11eb-3e6c-1324a5470226
# ╠═7b878c70-5234-11eb-39dd-9dc1bbfd76f3
# ╟─7a5722c0-5234-11eb-21c9-ef24bfde97b2
# ╠═7c1cc790-5234-11eb-0c3e-0381bf9a6b02
# ╟─989fab30-5234-11eb-1360-176b9497415b
# ╠═9758d350-5234-11eb-0470-7de7fb74ebbd
# ╟─f1371bb0-523f-11eb-37cb-c1b9ab9ee329
# ╠═218e2380-5240-11eb-27e7-374a37c7100f
# ╟─43c4d610-5240-11eb-3bbf-993a98d50ea5
# ╠═4dcac110-5240-11eb-2ce5-df407e48f332
# ╟─b7c48182-5242-11eb-2eea-7f7090854eaa
# ╠═b31180c0-5242-11eb-21a0-1175781341f2
# ╟─7ac2b430-5294-11eb-180b-c7e8a015731a
# ╟─9aa625f0-52aa-11eb-0ed6-eb5a4723e806
# ╟─97b5a370-52af-11eb-2029-ff804f96a261
# ╠═103a8d10-52ab-11eb-14a5-8915384a290b
# ╠═5d1c19a0-52ab-11eb-044a-2f0d527709f5
# ╟─809e81e0-52b2-11eb-22f5-e12d7e2a10b2
# ╟─5c0d2430-52b7-11eb-1a30-ed8986c0a8fb
# ╠═90734ba0-52b7-11eb-0723-652371f58cc1
# ╠═a6e01c10-52b7-11eb-1366-05747462a345
# ╠═ab5de9c2-52bc-11eb-1186-292123827b72
# ╠═abf1772e-52bc-11eb-33b4-c39608653da6
# ╟─03a3f4fe-52bf-11eb-3fe5-1d21d702651d
# ╠═55702992-52be-11eb-07b9-b57e098f7b3c
# ╟─733df1d0-52c0-11eb-2217-43b21eed1187
# ╠═d6f36a82-52bf-11eb-294b-e5276cabc27b
# ╟─2bf99490-52e4-11eb-0ddb-2b9c799865e7
# ╟─4b66edf0-52e4-11eb-25ab-0978549525c6
