abstract EnKF_iter 
typealias IntPairs Tuple{Int64, Int64}

immutable EAKF <: EnKF_iter end


function update( iiz :: Array{IntPairs}, loc :: Localization, kf :: EnKF_iter, 
                 rraw :: Array{Float64,2},
                 obs:: Array{Float64,1},error_obs :: Array{Float64, 1})
  res = copy(rraw) 
  for ii = iiz
    (ob, err_ob, ra, y) = iter_select(ii, loc, res, obs, error_obs)
    upyy = iter_update(ii, loc, kf, ra, y, ob, err_ob)
    iter_refill!(ii, loc, res, upyy)
  end
  return res
end

function iter_update( ii :: IntPairs, loc :: NearPoints, kf :: EAKF, 
                      raw0 :: Array{Float64,2}, y0 :: Array{Float64,1}, 
                      obs :: Float64, error_obs :: Float64, 
                    )
  y1 = y0 - mean(y0)
  raw1 = raw0 .- mean(raw0,2)
  covr = (raw1 * y1) / sum(y1.^2)

  ## Updating the mean
  sigma_ratio = err_obs^2 / (error_obs^2 + var(y0))
  delta_ym = (obs - mean(y0) ) * (1-sigma_ratio)
  delta_yd :: Array{Float64,1} = y1 (sqrt(sigma_ratio)-1)
  delta_upyy0 :: Array{Float64,1} = covr* delta_y
  delta_upyy1 :: Array{Float64,2} = covr * delta_yd'
  return ((delta_upyy0 .+ delta_upyy1) + raw0)
end

