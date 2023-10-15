\name{ref}
\alias{ref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a reference to an object
}
\description{
Create a reference to an object. Nature of this reference is an on-the-fly function. Use \code{rm()} to remove it.
}
\usage{
ref(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A character string as the alias.
}
  \item{y}{
A character string, name of the object to be referenced.
}
}
\details{
This function makes use of \code{makeActiveBinding()}. It will search for the name \code{y} from the current scope (environment) and up. Once found, \code{x} will work as a reference to \code{y}.

If \code{y} is something like \code{"tmp$c"}, and let us assume there is a \code{tmp = list(a = 1:2, b = 3:4)} in the current scope, and also a \code{tmp = list(c = 5:6, d = 7:8)} in the parent scope, then the function will stop at the current environment and throw an error, because \code{tmp} in the current scope has no member \code{$c}.
}
\value{
None
}
\examples{
# tmp = list(a = 1:3, b = 2:4, c = 2:5)
# 
# f = function()
# {
#   tmp = list(a = 11:12, b = 13:15)
#   ref("x", "tmp$c")
#   x = 1:10
#   print(tmp)
#   # rm(x) or deref("x"), x will not go outside function scope.
# }
# 
# f(); tmp; x # f() error; tmp unchanged; x does not exist.
# 
# g = function()
# {
#   tmp = list(a = 11:12, b = 13:15)
#   ref("x", "tmp$a")
#   x = 1:10
#   print(tmp)
# }
# g()
}
















