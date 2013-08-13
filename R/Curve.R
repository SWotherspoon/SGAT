curve.model <- function(datetime, light, segments,
                        twilight.model = c("LogNormal"),
                        alpha, beta,
                        logp.x = function(x) rep.int(0L, nrow(x)),
                        logp.z = function(z) rep.int(0L, nrow(z)),
                        x0, z0 = NULL,
                        fixedx = FALSE) {




}

