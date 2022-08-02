#' ---
#' title: Processamento de um inventário florestal utilizando funções de afilamento  
#' subtitle: Trabalho 2
#' author: Gabriel F. Pereira
#' ---

arv <- read.csv2( "exercícios/exercício_02/cubagem.csv" )
arv$dap <- arv$cap / pi
arv$disc <- arv$dicc - ( 2 * arv$espcasca )
espcasca_dap <- subset( arv, hi == 1.3 )
espcasca_dap$dapsc <- espcasca_dap$dap - ( 2 * espcasca_dap$espcasca )
espcasca_dap <- espcasca_dap[ c( "arvore", "dapsc" ) ]
arv <- merge( arv, espcasca_dap, by = "arvore" )

par( mfrow = c( 1, 2 ),
     pch = 20,
     bg = "#171717",
     col.axis = "white",
     col.lab = "white",
     col.main = "white",
     fg = "white" ) 
with( arv, plot( hi/ht, dicc, pch = 20 ) )
with( arv, plot( hi, dicc, pch = 20 ) )


#' # Ajuste 5º grau ----------------------------------------------------------------

#' #### Ajuste com casca

#' ##### Ajuste linear
ajlin <- lm( formula = "I(dicc/dap)~I(hi/ht)+I((hi/ht)^2)+I((hi/ht)^3)+I((hi/ht)^4)+I((hi/ht)^5)", data = arv )
summary( ajlin )
bs <- coef( ajlin )

#' ##### Ajuste não-linear
ajnlin <- nls( formula = "dicc~dap*(b0+b1*(hi/ht)+b2*((hi/ht)^2)+b3*((hi/ht)^3)+b4*((hi/ht)^4)+b5*((hi/ht)^5))",
               data = arv,
               start = list( b0 = bs[1], b1 = bs[2], b2 = bs[3], b3 = bs[4], b4 = bs[5], b5 = bs[6] ) )

#' ##### Erro padrão residual
syx <- summary( ajnlin )$sigma
syxperc <- syx/mean( arv$dicc ) * 100


#' #### Ajuste sem casca 

#' ##### Ajuste linear
ajlinsc <- lm(formula = "I(disc/dapsc)~I(hi/ht)+I((hi/ht)^2)+I((hi/ht)^3)+I((hi/ht)^4)+I((hi/ht)^5)", data = arv )
summary( ajlinsc )
bs <- coef( ajlinsc )

#' ##### Ajuste não linear
ajnlinsc <- nls( "disc~dapsc*(b0+b1*(hi/ht)+b2*((hi/ht)^2)+b3*((hi/ht)^3)+b4*((hi/ht)^4)+b5*((hi/ht)^5))",
                 data = arv,
                 start = list( b0 = bs[1], b1 = bs[2], b2 = bs[3], b3 = bs[4], b4 = bs[5], b5 = bs[6] ) )

#' ##### Erro padrão residual
syxsc <- summary( ajnlinsc )$sigma
syxpercsc <- syxsc/mean( arv$disc ) * 100


#' # Ajuste Methol  ----------------------------------------------------------------

#' #### Ajuste com casca

#' ##### Ajuste não linear
methol <- nls( "dicc~(a0*(dap^a1)*(a2^dap))*(1-(sqrt(hi/ht)))^(b1*log(hi/ht+0.001)+(b2*exp(hi/ht))+(b3*(dap/ht))+(b4*log(dap))+(b5*(ht/sqrt(hi)))+(b6*((dap/ht)/hi)))",
               arv,
               start = list( a0 = 1.28423,
                             a1 = 1.05171,
                             a2 = 0.99141,
                             b1 = -0.28784,
                             b2 = 0.22426,
                             b3 = 0.12551,
                             b4 = -0.00775,
                             b5 = -0.006264,
                             b6 = -0.15257 ) )

#' ##### Erro padrão residual
mtsyx <- summary( methol )$sigma
mtsyxperc <- mtsyx/mean( arv$dicc ) * 100

#' #### Ajuste sem casca 

#' ##### Ajuste não linear
metholsc <- nls( "disc~(a0*(dapsc^a1)*(a2^dapsc))*(1-(sqrt(hi/ht)))^(b1*log(hi/ht+0.001)+(b2*exp(hi/ht))+(b3*(dapsc/ht))+(b4*log(dapsc))+(b5*(ht/sqrt(hi)))+(b6*((dapsc/ht)/hi)))",
                 arv,
                 start = list( a0 = 1.28423,
                               a1 = 1.05171,
                               a2 = 0.99141,
                               b1 = -0.28784,
                               b2 = 0.22426,
                               b3 = 0.12551,
                               b4 = -0.00775,
                               b5 = -0.006264,
                               b6 = -0.15257 ) )

#' ##### Erro padrão residual
mtsyxsc <- summary( metholsc )$sigma
mtsyxpercsc <- mtsyxsc/mean( arv$disc ) * 100


#' # Ajuste Kozak  ----------------------------------------------------------------

#' #### Ajuste com casca

#' ##### Ajuste não linear
kozak <- nls( "dicc~(a0*(dap^a1)*(ht^a2))*((1-(hi/ht)^(1/3))/(1-(1.3/ht)^(1/3)))^((b1*(hi/ht)^4)+(b2*(1/(exp(dap/ht))))+(b3*((1-(hi/ht)^(1/3))/(1-(1.3/ht)^(1/3)))^0.1)+(b4*(1/dap))+(b5*ht^(1-(hi/ht)^1/3))+(b6*((1-(hi/ht)^(1/3))/(1-(1.3/ht)^(1/3)))))",
              arv,
              start = list(a0 = 1.17996,
                           a1 = 0.93279,
                           a2 = 0.01165,
                           b1 = 0.25256,
                           b2 = -0.62129,
                           b3 = 0.711660,
                           b4 = 2.075633,
                           b5 = 0.040055,
                           b6 = -0.31407 ) )

#' ##### Erro padrão residual
kosyx <- summary( kozak )$sigma
kosyxperc <- kosyx/mean( arv$dicc ) * 100

#' #### Ajuste sem casca

#' ##### Ajuste não linear
kozaksc <- nls( "disc~(a0*(dapsc^a1)*(ht^a2))*((1-(hi/ht)^(1/3))/(1-(1.3/ht)^(1/3)))^((b1*(hi/ht)^4)+(b2*(1/(exp(dapsc/ht))))+(b3*((1-(hi/ht)^(1/3))/(1-(1.3/ht)^(1/3)))^0.1)+(b4*(1/dapsc))+(b5*ht^(1-(hi/ht)^1/3))+(b6*((1-(hi/ht)^(1/3))/(1-(1.3/ht)^(1/3)))))",
                arv,
                start = list(a0 = 1.17996,
                             a1 = 0.93279,
                             a2 = 0.01165,
                             b1 = 0.25256,
                             b2 = -0.62129,
                             b3 = 0.711660,
                             b4 = 2.075633,
                             b5 = 0.040055,
                             b6 = -0.31407 ) )

#' ##### Erro padrão residual
kosyxsc <- summary( kozaksc )$sigma
kosyxpercsc <- kosyxsc/mean( arv$disc ) * 100


#' # Ajuste Ormerod  ----------------------------------------------------------------

#' #### Ajuste com casca

#' ##### Ajuste não linear
ormerod <- nls("dicc~dap*((ht-hi)/(ht-1.3))^b1",
               arv,
               start = list( b1 = 0.1 ) )

#' ##### Erro padrão residual
orsyx <- summary( ormerod )$sigma
orsyxperc <- orsyx/mean( arv$dicc ) * 100

#' #### Ajuste sem casca
ormerodsc <- nls("disc~dapsc*((ht-hi)/(ht-1.3))^b1",
                 arv,
                 start = list( b1 = 0.1 ) )

#' ##### Erro padrão residual
orsyxsc <- summary( ormerodsc )$sigma
orsyxpercsc <- orsyxsc/mean( arv$disc ) * 100

#' # Comparação erro padrão residual

setNames( data.frame( Modelos = c( "Schöepfer (com casca)", "Schöepfer (sem casca)",
                                   "Methol (com casca)", "Methol (sem casca)",
                                   "Kozak 2004 (com casca)", "Kozak 2004 (sem casca)",
                                   "Ormerod (com casca)", "Ormerod (sem casca)" ), 
                      Erro = c( syxperc, syxpercsc, 
                                mtsyxperc, mtsyxpercsc,
                                kosyxperc, kosyxpercsc,
                                orsyxperc, orsyxpercsc)
                      
), c( "|Modelos|", "|Erro Padrão Residual|" ) )


#' #### Valores preditos

#' ##### 5° grau
arv$diccest5grau <- predict( ajnlin )
arv$discest5grau <- predict( ajnlinsc )

#' ##### methol
arv$diccestmet <- predict( methol )
arv$discestmet <- predict( metholsc )

#' ##### kozak
arv$diccestkoz <- predict( kozak )
arv$discestkoz <- predict( kozaksc )

#' ##### ormerod
arv$diccestorm <- predict( ormerod )
arv$discestorm <- predict( ormerodsc )


#' #### Resíduos

#' ##### 5° grau
arv$res5grau <- residuals( ajnlin )
arv$res5grausc <- residuals( ajnlinsc )

#' ##### methol
arv$resmet <- residuals( methol )
arv$resmetsc <- residuals( metholsc )

#' ##### kozak
arv$reskoz <- residuals( kozak )
arv$reskozsc <- residuals( kozaksc )

#' ##### ormerod
arv$resorm <- residuals( ormerod )
arv$resormsc <- residuals( ormerodsc )


#' # Análise Gráfica   ----------------------------------------------------------------

library( fBasics )

#' ### Com casca
par( mfrow = c( 2, 2 ),
     pch = 20,
     bg = "#171717",
     col.axis = "white",
     col.lab = "white",
     col.main = "white",
     fg = "white" ) 

qqnormPlot( ( arv$res5grausc ), title = FALSE, main = "Schöepfer" )
qqnormPlot( ( arv$resmetsc ), title = FALSE, main = "Methol" )
qqnormPlot( ( arv$reskozsc ), title = FALSE, main = "Kozak" )
qqnormPlot( ( arv$resormsc ), title = FALSE, main = "Ormerod" )

graficos_res <- function( hi, res, modelo, xlab = "hi(m)" ) {
  plot( hi, res, 
        pch = 20, 
        ylim = c( -5, 5 ),
        ylab = "Resíduo (m³)", 
        main = modelo )
  abline( h = 0, col = "red" )
}

graficos_res( arv$hi, arv$res5grau, "Schöepfer" )
graficos_res( arv$hi, arv$resmet, "Methol" )
graficos_res( arv$hi, arv$reskoz, "Kozak" )
graficos_res( arv$hi, arv$resorm, "Ormerod" )

#' ### Sem casca
par( mfrow = c( 2, 2 ),
pch = 20,
bg = "#171717",
col.axis = "white",
col.lab = "white",
col.main = "white",
fg = "white" ) 

qqnormPlot((arv$res5grausc),title = FALSE, main = "Schöepfer")
qqnormPlot((arv$resmetsc),title = FALSE, main = "Methol")
qqnormPlot((arv$reskozsc),title = FALSE, main = "Kozak")
qqnormPlot((arv$resormsc),title = FALSE, main = "Ormerod")

graficos_res( arv$hi, arv$res5grausc, "Schöepfer" )
graficos_res( arv$hi, arv$resmetsc, "Methol" )
graficos_res( arv$hi, arv$reskozsc, "Kozak" )
graficos_res( arv$hi, arv$resormsc, "Ormerod" )
