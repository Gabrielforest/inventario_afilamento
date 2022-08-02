#' ---
#' title: Função multiprodutos 
#' subtitle: Alterando para outros modelos
#' author: Gabriel F. Pereira
#' ---

library( cmrinvflor )
??cmrinvflor::multprodarvbt5grau

produtos <- openxlsx::read.xlsx( "aulas/aula_02/aplic_5grau_r.xlsx", sheet = "produtos" )

#' ## Parâmetros da função
fustes <- openxlsx::read.xlsx( "aulas/aula_02/aplic_5grau_r.xlsx", sheet = "fustes" )
coefs <- openxlsx::read.xlsx( "aulas/aula_02/aplic_5grau_r.xlsx", sheet = "coefs" )
nomprod <- as.matrix( produtos[, 1] )
vprod <- as.matrix( produtos[, 2:5] )
subseq <- 1

procafill <- as.data.frame( cmrinvflor::multprodarvbt5grau( fustes, coefs, vprod, nomprod, subseq ) )

#' #### Dados importantes que são gerados
lapply( procafill[, 23:28 ], sum )

d5grau <- function( hi, dap, htafil, 
                    b0, b1, b2, b3, b4, b5 ) {
  if( hi >= htafil ) { 0 
  } else { dap * (b0 + b1 * hi/ htafil + b2 * (( hi/ htafil )^2) 
                  + b3 * ((hi / htafil )^3 ) + b4 * ((hi/htafil)^4) + 
                    b5 * ((hi/ htafil)^5) ) }
}

#' #### Alterando para multiprodutos usando Kozak

#' Substituindo função v5grau pela função vkozak
v5grau <- function( htoco, dap, b0, b1, b2, b3, b4, b5 ) {
  ( pi / 40000 ) * dap^2 * 
    ( b0^2 * htoco + b0 * b1 * htoco^2 + 
        ( 2/3 * b0 * b2 + 1/3 * b1^2 ) * htoco^3 + 
        ( 1/2 * b0 * b3 + 1/2 * b1 * b2 ) * htoco^4 +
        ( 2/5 * b0 * b4 + 2/5 * b1 * b3 + 1/5 * b2^2 ) * htoco^5 +
        ( 1/3 * b0 * b5 + 1/3 * b1 * b4 + 1/3 * b2 * b3) * htoco^6 +
        ( 2/7 * b1 * b5 + 2/7 * b2 * b4 + 1/7 * b3^2 ) * htoco^7 +
        ( 1/4 * b2 * b5 + 1/4 * b3 * b4) * htoco^8 +
        ( 2/9 * b3 * b5 + 1/9 * b4^2 ) * htoco^9 +
        1/5 * b4 * b5 * htoco^10 + 1/11 * b5^2 * htoco^11 )
}
v5grau( htoco, dap, b0, c1, c2, c3, c4, c5 )

#' ### Integral de Kozak
fgkozak <- function( hi, dap, ht, b0, b1, b2, b3, b4, b5, b6, b7, b8 ) {
  return(pi/40000*(b0*(dap**b1)*(ht**b2)*
                     (((1-((hi/ht)^(1/3)))/(1-((1.3/ht)^(1/3))))**(b3*(hi/ht)**4+
                                                                     b4*(1/exp(dap/ht))+(b5*(((1-((hi/ht)^(1/3)))/(1-((1.3/ht)^(1/3))))**0.1))+(b6*(1/dap))+
                                                                     (b7*(ht**(1-((hi/ht)**(1/3)))))+(b8*((1-((hi/ht)^(1/3)))/(1-((1.3/ht)^(1/3))))))))^2)
}
vkozak <- function( hi, dap, ht, b0, b1, b2, b3, b4, b5, b6, b7, b8 ){ 
  return(integrate(fgkozak,lower=0,upper=hi,dap=dap,ht=ht,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,b5=b5,b6=b6,b7=b7,b8=b8)$value) 
}


function (fustes, coefs, vprod, nomprod, subseq) 
{
  # dados são "corrigidos" e ajustados de acordo com medições
  fustes$htafil <- fustes$htre
  ii <- (fustes$haprov_total == 0 & fustes$htre < fustes$htest)
  fustes$htafil[ii] <- fustes$htest[ii]
  ii <- fustes$htafil < 1.3
  if (sum(ii) > 0) {
    fustes <- fustes[ii == FALSE, ]
  }
  if (is.vector(vprod)) 
    vprod <- t(vprod)
  nprod <- nrow(vprod)
  if (subseq == 1 & nprod > 1) {
    ii <- order(vprod[, 1], decreasing = TRUE)
    vprod <- vprod[ii, ]
    nomprod <- nomprod[ii]
  }
  # cria a tabela principal (dados) e as auxiliares que entrarão na principal
  nl <- nrow(fustes)
  ncols <- ncol(fustes)
  dados <- merge(fustes, coefs, by.x = "idestcoef", by.y = "idestcoef", 
                 sort = FALSE)
  rm(fustes, coefs)
  vnt <- matrix(NA, nl, 5)
  vest <- matrix(NA, nl, 1)
  vhd <- matrix(NA, nl, 1)
  vtoco <- matrix(NA, nl, 1)
  vntoras <- matrix(NA, nl, nprod)
  vntoras[, ] <- 0
  vsdtoras <- matrix(NA, nl, nprod)
  vsdtoras[, ] <- 0
  vhif <- matrix(NA, nl, nprod)
  vhif[, ] <- 0
  vdmin <- matrix(NA, nl, nprod)
  vdmin[, ] <- 0
  vvoltoras <- matrix(NA, nl, nprod + 1)
  vvoltoras[, ] <- 0
  # Contas começam aqui (estimativas usando polinômio do quinto grau)
  # alterar os betas
   dados$c0 <- dados$b0
   dados$c1 <- with(dados, b1/htafil)
   dados$c2 <- with(dados, b2/(htafil^2))
   dados$c3 <- with(dados, b3/(htafil^3))
   dados$c4 <- with(dados, b4/(htafil^4))
   dados$c5 <- with(dados, b5/(htafil^5))
  # não usamos subseq == 0
  if (subseq == 0) {
    vtocos <- matrix(NA, nl, nprod)
    vtocos[, ] <- 0
    htoco <- vprod[1, 4]
    vtoco[, 1] <- with(dados, v5grau(vprod[1, 4], dap, c0, 
                                     c1, c2, c3, c4, c5))
    for (np in 1:nprod) {
      vnt[, ] <- 0
      vest[, 1] <- 0
      vhd[, 1] <- vprod[np, 4]
      if (htoco != vprod[np, 4]) {
        htoco = vprod[np, 4]
        vtoco[, 1] <- with(dados, v5grau(htoco, dap, 
                                         c0, c1, c2, c3, c4, c5))
      }
      limite <- vprod[np, 1]
      aux <- 0
      while (aux == 0) {
        aux <- 1
        for (i in 1:nl) {
          if (vnt[i, 2] == 0) {
            if (dados$htre[i] >= (htoco + vprod[np, 2])) {
              if ((vhd[i, 1] + vprod[np, 2]) <= dados$htre[i]) {
                hc = vhd[i, 1] + (vprod[np, 2]/2)
                vhd[i, 1] = vhd[i, 1] + vprod[np, 2]
                di = with(dados[i, ], d5grau(vhd[i, 1], 
                                             dap, htafil, b0, b1, b2, b3, b4, b5))
                if (di >= limite) {
                  vnt[i, 1] = vnt[i, 1] + 1
                  vnt[i, 4] = di
                  vnt[i, 5] = vnt[i, 5] + with(dados[i, 
                  ], d5grau(hc, dap, htafil, b0, b1, 
                            b2, b3, b4, b5))
                  aux = 0
                }
                else {
                  vnt[i, 2] = 1
                  if (vprod[np, 3] == vprod[np, 2]) {
                    vnt[i, 3] = vprod[np, 2]
                  }
                  else {
                    inc <- vprod[np, 3]
                    vhd[i, 1] = vhd[i, 1] - vprod[np, 
                                                  2]
                    hc = vhd[i, 1]
                    hi = vhd[i, 1] + inc
                    if (hi <= dados$htre[i]) {
                      di = limite
                      while (di >= limite & hi <= dados$htre[i]) {
                        di = with(dados[i, ], d5grau(hi, 
                                                     dap, htafil, b0, b1, b2, b3, 
                                                     b4, b5))
                        if (di >= limite) {
                          inc <- 0.3
                          hi = hi + inc
                          vnt[i, 4] = di
                        }
                      }
                      hi = hi - inc
                      if (inc == vprod[np, 3]) {
                        vnt[i, 3] = vprod[np, 2]
                      }
                      else {
                        vnt[i, 1] = vnt[i, 1] + 1
                        vnt[i, 3] = hi - vhd[i, 1]
                        hc = hc + (vnt[i, 3]/2)
                        vnt[i, 5] = vnt[i, 5] + with(dados[i, 
                        ], d5grau(hc, dap, htafil, 
                                  b0, b1, b2, b3, b4, b5))
                      }
                    }
                    else {
                      if (vn[i, 1] > 0) {
                        vnt[i, 3] = vprod[np, 2]
                      }
                      else {
                        vnt[i, 3] = 0
                      }
                    }
                  }
                }
              }
              else {
                vnt[i, 2] = 1
                inc <- vprod[np, 3]
                hc = vhd[i, 1]
                hi = vhd[i, 1] + inc
                if (hi <= dados$htre[i]) {
                  di = limite
                  while (di >= limite & hi <= dados$htre[i] & 
                         di != 0) {
                    di = with(dados[i, ], d5grau(hi, 
                                                 dap, htafil, b0, b1, b2, b3, b4, 
                                                 b5))
                    if (di >= limite) {
                      inc <- 0.3
                      hi = hi + inc
                      vnt[i, 4] = di
                    }
                  }
                  hi = hi - inc
                  if (inc == vprod[np, 3]) {
                    vnt[i, 3] = vprod[np, 2]
                  }
                  else {
                    vnt[i, 1] = vnt[i, 1] + 1
                    vnt[i, 3] = hi - vhd[i, 1]
                    hc = hc + (vnt[i, 3]/2)
                    vnt[i, 5] = vnt[i, 5] + with(dados[i, 
                    ], d5grau(hc, dap, htafil, b0, 
                              b1, b2, b3, b4, b5))
                  }
                }
                else {
                  vnt[i, 3] = vprod[np, 2]
                }
              }
            }
            else {
              vnt[i, 2] <- 1
              if (dados$htre[i] >= (htoco + vprod[np, 
                                                  3])) {
                hc <- vhd[i, 1]
                hi <- vhd[i, 1] + vprod[np, 3]
                di <- limite
                while ((di >= limite) & (hi < dados$htre[i])) {
                  di <- with(dados[i, ], d5grau(hi, dap, 
                                                htafil, b0, b1, b2, b3, b4, b5))
                  if (di >= limite) {
                    hi <- hi + 0.3
                    vnt[i, 4] <- di
                  }
                }
                hi <- hi - (vhd[i, 1] + 0.3)
                if (hi < vprod[np, 3]) {
                  vnt[i, 3] = vprod[np, 2]
                }
                else {
                  vnt[i, 1] = vnt[i, 1] + 1
                  vnt[i, 3] = hi
                  hc = hc + (vnt[i, 3]/2)
                  vnt[i, 5] = vnt[i, 5] + with(dados[i, 
                  ], d5grau(hc, dap, htafil, b0, b1, 
                            b2, b3, b4, b5))
                }
              }
            }
          }
        }
      }
      for (i in 1:nl) {
        vtocos[i, np] = vtoco[i]
        if (vnt[i, 1] > 0) {
          vntoras[i, np] = vnt[i, 1]
          vsdtoras[i, np] = vnt[i, 5]
          vhif[i, np] = vnt[i, 3]
          vdmin[i, np] = vnt[i, 4]
          vhd[i, 1] = vprod[np, 4] + ((vnt[i, 1] * vprod[np, 
                                                         2]) - (vprod[np, 2] - vnt[i, 3]))
          vest[i, 1] <- with(dados[i, ], v5grau(vhd[i, 
                                                    1], dap, c0, c1, c2, c3, c4, c5))
          if (vest[i, 1] > 0) {
            vvoltoras[i, np] = vest[i, 1] - vtoco[i]
          }
        }
      }
    }
    vvoltoras[, nprod + 1] <- with(dados, v5grau(htre, dap, 
                                                 c0, c1, c2, c3, c4, c5))
    calcarv <- as.matrix(cbind(dados$idfustemed, vntoras, 
                               vsdtoras, vdmin, vhif, vtocos, vvoltoras))
    nome <- c("idfustemed")
    for (i in 1:nprod) {
      nome <- c(nome, paste("ntorasi", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("sdtorasi", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("dmini", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("toramini", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("vtocoi", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("vprodi", nomprod[i], 
                            sep = ""))
    }
    nome <- c(nome, "vti")
    dimnames(calcarv) <- list(c(), nome)
  }
  # subseq == 1
  else {
    vhvolac <- matrix(NA, nl, 2)
    vhvolac[, ] <- 0
    vhd[, 1] <- vprod[1, 4]
    vhvolac[, 1] <- vprod[1, 4]
    # alterar v5grau 
    vhvolac[, 2] <- with(dados, v5grau(vhvolac[, 1], dap, 
                                       c0, c1, c2, c3, c4, c5))
    vtoco[, 1] <- vhvolac[, 2]
    for (np in 1:nprod) {
      vnt[, ] <- 0
      vest[, 1] <- 0
      limite <- vprod[np, 1]
      aux = 0
      while (aux == 0) {
        aux <- 1
        for (i in 1:nl) {
          if (vnt[i, 2] == 0) {
            if ((vhd[i, 1] + vprod[np, 2]) <= dados$htre[i]) {
              hc = vhd[i, 1] + (vprod[np, 2]/2)
              vhd[i, 1] = vhd[i, 1] + vprod[np, 2]
              di = with(dados[i, ], d5grau(vhd[i, 1], 
                                           dap, htafil, b0, b1, b2, b3, b4, b5))
              # parou aq
              if (di >= limite) {
                vnt[i, 1] = vnt[i, 1] + 1
                vnt[i, 4] = di
                vnt[i, 5] = vnt[i, 5] + with(dados[i, 
                ], d5grau(hc, dap, htafil, b0, b1, 
                          b2, b3, b4, b5))
                aux = 0
              }
              else {
                vnt[i, 2] = 1
                if (vprod[np, 3] == vprod[np, 2]) {
                  vnt[i, 3] = vprod[np, 2]
                }
                else {
                  inc <- vprod[np, 3]
                  vhd[i, 1] = vhd[i, 1] - vprod[np, 2]
                  hc = vhd[i, 1]
                  hi = vhd[i, 1] + inc
                  if (hi <= dados$htre[i]) {
                    di = limite
                    while (di >= limite & hi <= dados$htre[i]) {
                      di = with(dados[i, ], d5grau(hi, dap, htafil, 
                                                   b0, b1, b2, b3, b4, b5))
                      if (di >= limite) {
                        inc <- 0.3
                        hi = hi + inc
                        vnt[i, 4] = di
                      }
                    }
                    hi = hi - inc
                    if (inc == vprod[np, 3]) {
                      vnt[i, 3] = vprod[np, 2]
                    }
                    else {
                      vnt[i, 1] = vnt[i, 1] + 1
                      vnt[i, 3] = hi - vhd[i, 1]
                      hc = hc + (vnt[i, 3]/2)
                      vnt[i, 5] = vnt[i, 5] + with(dados[i, 
                      ], d5grau(hc, dap, htafil, 
                                b0, b1, b2, b3, b4, b5))
                    }
                  }
                  else {
                    vnt[i, 3] = vprod[np, 2]
                  }
                }
              }
            }
            else {
              vnt[i, 2] = 1
              inc <- vprod[np, 3]
              hc = vhd[i, 1]
              hi = vhd[i, 1] + inc
              if (hi <= dados$htre[i]) {
                di = limite
                while (di >= limite & hi <= dados$htre[i] & 
                       di != 0) {
                  di = with(dados[i, ], d5grau(hi, dap, htafil, 
                                               b0, b1, b2, b3, b4, b5))
                  if (di >= limite) {
                    inc <- 0.3
                    hi = hi + inc
                    vnt[i, 4] = di
                  }
                }
                hi = hi - inc
                if (inc == vprod[np, 3]) {
                  vnt[i, 3] = vprod[np, 2]
                }
                else {
                  vnt[i, 1] = vnt[i, 1] + 1
                  vnt[i, 3] = hi - vhd[i, 1]
                  hc = hc + (vnt[i, 3]/2)
                  vnt[i, 5] = vnt[i, 5] + with(dados[i, 
                  ], d5grau(hc, dap, htafil, b0, b1, 
                            b2, b3, b4, b5))
                }
              }
              else {
                vnt[i, 3] = vprod[np, 2]
              }
            }
          }
        }
      }
      for (i in 1:nl) {
        if (vnt[i, 1] > 0) {
          vntoras[i, np] = vnt[i, 1]
          vsdtoras[i, np] = vnt[i, 5]
          vhif[i, np] = vnt[i, 3]
          vdmin[i, np] = vnt[i, 4]
          hi = (vnt[i, 1] * vprod[np, 2]) - (vprod[np, 
                                                   2] - vnt[i, 3])
          vhvolac[i, 1] = vhvolac[i, 1] + hi
          vest[i, 1] <- with(dados[i, ], v5grau(vhvolac[i, 
                                                        1], dap, c0, c1, c2, c3, c4, c5))
          if (vest[i, 1] > 0) {
            vvoltoras[i, np] = vest[i, 1] - vhvolac[i, 
                                                    2]
            vhvolac[i, 2] = vest[i, 1]
          }
        }
        vhd[i, 1] = vhvolac[i, 1]
      }
    }
    vvoltoras[, nprod + 1] <- with(dados, v5grau(htre, dap, 
                                                 c0, c1, c2, c3, c4, c5))
    calcarv <- as.matrix(cbind(dados$idfustemed, vntoras, 
                               vsdtoras, vdmin, vhif, vtoco[, 1], vvoltoras))
    nome <- c("idfustemed")
    for (i in 1:nprod) {
      nome <- c(nome, paste("ntorasi", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("sdtorasi", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("dmini", nomprod[i], 
                            sep = ""))
    }
    for (i in 1:nprod) {
      nome <- c(nome, paste("toramini", nomprod[i], 
                            sep = ""))
    }
    nome <- c(nome, "vtocoi")
    for (i in 1:nprod) {
      nome <- c(nome, paste("vprodi", nomprod[i], 
                            sep = ""))
    }
    nome <- c(nome, "vti")
    dimnames(calcarv) <- list(c(), nome)
  }
  return(calcarv)
}