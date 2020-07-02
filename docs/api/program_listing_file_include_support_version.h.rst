
.. _program_listing_file_include_support_version.h:

Program Listing for File version.h
==================================

|exhale_lsh| :ref:`Return to documentation for file <file_include_support_version.h>` (``include/support/version.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   /* #############################################
    *             This file is part of
    *                    ZERO
    *
    *             Copyright (c) 2020
    *     Released under the Creative Commons
    *        Zero v1.0 Universal License
    *
    *              Find out more at
    *        https://github.com/ds4dm/ZERO
    * #############################################*/
   
   
   #include <string>
   #define ZERO_VERSION "1.0.0"
   #define ZERO_VERSION_MAJOR "1"
   #define ZERO_VERSION_MINOR "0"
   #define ZERO_VERSION_PATCH "0"
   
   class ZEROVersion {
   public:
     ZEROVersion(std::string &major, std::string &minor, std::string &patch) {
        major = ZERO_VERSION_MAJOR;
        minor = ZERO_VERSION_MINOR;
        patch = ZERO_VERSION_PATCH;
     }
     ZEROVersion(std::string &version) { version = ZERO_VERSION; }
   };
