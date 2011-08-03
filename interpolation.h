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
                // n - ���������� ������� �������� � ����� ������ ������ 3
                Matrix T(unsigned int n,double t);
                Matrix V(unsigned int n,double k);
                // N - ���������� ���������
                //
                Matrix Q(unsigned int N);
                // �����
//                std::vector<Matrix> U;
        public:
                void initializeWith(std::vector<Matrix> U);
                double calculate(double t);
};

}
//---------------------------------------------------------------------------
#endif
