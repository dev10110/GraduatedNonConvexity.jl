module GraduatedNonConvexity

using Printf

include("utils.jl")
include("gnc_gm.jl")
include("gnc_tls.jl")

export GNC_GM, GNC_TLS
export GNC_GM!, GNC_TLS!

end
