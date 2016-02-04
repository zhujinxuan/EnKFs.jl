abstract Localization

include("dim_select.jl")
export iter_dim_seletct

function iter_select (ii :: Tuple{Int64, Int64}, loc :: Localization,
                      raw :: Array{Float64,2},
                      obs :: Array{Float64,1}, error_obs :: Array{Float64,1})
  (iobs, ixx ) = ii
  llrr = iter_dim_seletct(ii, loc, size(raw))
  return (obs[iobs], error_obs[iobs], raw[llrr,:], squeeze(raw[ixx,:],1))
end

function iter_refill!(ii :: Tuple{Int64, Int64}, m :: Localization,
                     raw :: Array{Float64,2}, upyy :: Array{Float64,2} )
  (iobs, ixx ) = ii
  llrr = iter_dim_seletct(ii, loc, size(raw))
  raw[llrr,:] = upyy
end


export iter_select, iter_refill!
