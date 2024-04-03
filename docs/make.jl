using Wquantiles
import Pkg

Pkg.add("Documenter")
using Documenter

makedocs(
	sitename = "Wquantiles",
	format = Documenter.HTML(),
	modules = [Wquantiles]
	)

	# Documenter can also automatically deploy documentation to gh-pages.
	# See "Hosting Documentation" and deploydocs() in the Documenter manual
	# for more information.
	deploydocs(
		repo = "github.com/scottrsm/Wquantiles.jl.git"
	)
