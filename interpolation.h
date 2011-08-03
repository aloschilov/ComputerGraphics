//---------------------------------------------------------------------------

#ifndef interpolationH
#define interpolationH
#include "matrix.h"
#include <vector>

namespace Aloschil
{

class Spline
{
        protected:
                // n - определяет степень полинома в нашем случае всегда 3
                Matrix T(unsigned int n,double t);
                Matrix V(unsigned int n,double k);
                // N - количество сегментов
                //
                Matrix Q(unsigned int N);
                // точки
//                std::vector<Matrix> U;
        public:
                void initializeWith(std::vector<Matrix> U);
                double calculate(double t);
};

}
//---------------------------------------------------------------------------
#endif
