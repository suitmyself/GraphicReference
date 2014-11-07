/******************************************************************************
*   MyAssert.h
*
*   My own assert() macro which will give a message box but not terminate.
******************************************************************************/

#ifndef MYASSERT_H
#define MYASSERT_H

#include <stdio.h>
#include <windows.h>

#pragma option -wmsg
#ifndef NO_ASSERT
    #pragma message Asserts are ON. Define NO_ASSERT to disable
    #define assert(f)                                                                               \
        if (!(f)) {                                                                                 \
            char buf[80];                                                                           \
            sprintf(buf, "File : %s\nLine : %d", __FILE__, __LINE__);                               \
            MessageBox(0, buf, "Assert failed", MB_ICONWARNING);                                    \
        }
#else
    #pragma message Asserts are OFF. Undefine NO-ASSERT to enable
    #define assert(f)
#endif

#endif

