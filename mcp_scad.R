#### MCP and SCAD ####
sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))
  if(xn==0){
    return(x)
  }else{
    thval <- 1 - th/xn
    return(thval*((thval) >0)*x)
  }
}
mcp <- function(x,lam,nu,gam)
{
  temp <- gam*lam
  xn <- sqrt(sum(x^2))
  if(xn <= temp)
  {
    z <- sfun(x,lam/nu) / (1 - 1/(gam*nu))
  }else{
    z <- x
  }
  return(z)
}

scad <- function(x,lam,nu,gam)
{
  temp1 <- lam/nu
  temp2 <- gam * lam
  xn <- sqrt(sum(x^2))
  
  if(xn <= lam + temp1)
  {
    z <- sfun(x, temp1)
  }else if(xn <= temp2 & xn >= lam + temp1)
  {
    z <- sfun(x, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu))
  }else{
    z <- x
  }
  
  return(z)
  
}