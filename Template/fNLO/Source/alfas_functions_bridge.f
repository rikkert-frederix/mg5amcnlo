c Preserve the ALPHAS contract available to Python-generated model code.

      double precision function alphas(q)
      use alfas_functions_module, only: module_alphas => alphas
      implicit none
      double precision q

      alphas=module_alphas(q)
      end
