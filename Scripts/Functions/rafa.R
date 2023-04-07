###### wrapper para calcular func vetorizado

funcao_vetorizada <- function(x){
  result <- tryCatch(sua_func(x["a"],x["b"],x["c"],...), error = function(err) NA)
  return(result)
}
#### lista de parametros para funcao supostamente em dataframe e parametros na linha
resultado <- apply(df, 1, funcao_vetorizada, simplify = F)

resultado <- resultado[!is.na(resultado)] ### remove os NA

bind_resultado <- bind_rows(resultado, .id = "sample") ## aqui pra funcionar provavelmente precisa de elementos nomeados na lista



