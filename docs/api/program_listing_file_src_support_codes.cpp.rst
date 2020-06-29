
.. _program_listing_file_src_support_codes.cpp:

Program Listing for File codes.cpp
==================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_support_codes.cpp>` (``src/support/codes.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include "support/codes.h"
   #include <string>
   
   std::string std::to_string(const ZEROErrorCode &code) {
     string wrap = "ZEROErrorCode ";
     switch (code) {
     case ZEROErrorCode::MemoryError: {
       return wrap + "100: Memory error";
     }
     case ZEROErrorCode::InvalidQuery: {
       return wrap + "101: The attribute/data is not available";
     }
     case ZEROErrorCode::InvalidData: {
       return wrap + "102: The input data is invalid";
     }
     case ZEROErrorCode::SolverError: {
       return wrap + "103: A third-party solver has thrown an error";
     }
     case ZEROErrorCode::OutOfRange: {
       return wrap + "104: Index out of range";
     }
     case ZEROErrorCode::Numeric: {
       return wrap + "105: Numerical issues";
     }
     case ZEROErrorCode::IOError: {
       return wrap + "106: IO Interface error";
     }
     case ZEROErrorCode::Assertion: {
       return wrap + "107: Assertion failed";
     }
     default:
       return wrap + "0: Unknown error";
     }
   }
   std::string std::to_string(ZEROStatus st) {
     switch (st) {
     case ZEROStatus::NashEqNotFound:
       return std::string("NO_NASH_EQ_FOUND");
     case ZEROStatus::NashEqFound:
       return std::string("NASH_EQ_FOUND");
     case ZEROStatus::Solved:
       return std::string("SOLVED");
     case ZEROStatus::NotSolved:
       return std::string("NOT_SOLVED");
     case ZEROStatus::TimeLimit:
       return std::string("TIME_LIMIT");
     case ZEROStatus::Uninitialized:
       return std::string("UNINITIALIZED");
     case ZEROStatus::Numerical:
       return std::string("NUMERICAL_ISSUES");
     default:
       return std::string("UNKWNOWN");
     }
   }
