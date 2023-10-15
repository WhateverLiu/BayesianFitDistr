

# x and y are character strings.
ref <- function(x, y)
{

  tryCatch(eval(parse(text = paste0("tracemem(", y, ")")), envir = parent.frame()), error = function(e) { stop(paste0("keyALGs::ref() cannot find ", y)) })
  eval(parse(text = paste0("makeActiveBinding('", x, "', function(u) if(missing(u)) ", y, " else ", y, " <<- u, environment())")), envir = parent.frame())
}


deref <- function(x)
{
  rm(x, envir = parent.frame())
}


if(F)
{


  tmp = list(a = 1:3, b = 2:4, c = 2:5)

  f = function()
  {
    tmp = list(a = 11:12, b = 13:15)
    ref("x", "tmp$c")
    x = 1:10
    print(tmp)
    # rm(x) or deref("x"), x will not go outside function scope.
  }

  f(); tmp; x # f() error; tmp unchanged; x does not exist.

  g = function()
  {
    tmp = list(a = 11:12, b = 13:15)
    ref("x", "tmp$a")
    x = 1:10
    print(tmp)
  }
  g()


}










