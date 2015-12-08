#include <oxstd.h>

extern "test,FnFortranTest" FnFortranTest(const aflag, const adouble, const aint);

main()
{
    decl fl, i, d;

    FnFortranTest(&fl, &d, &i);

    print("\nboolean=", fl, " double=", d, " int=", i, "\n");
}

