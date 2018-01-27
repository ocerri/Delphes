#define private public
#define protected public
/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PileUpMergerPythia8.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PileUpMergerPythia8+;

#endif
// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME tmpdImodulesdIPythia8Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "modules/PileUpMergerPythia8.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PileUpMergerPythia8(void *p = 0);
   static void *newArray_PileUpMergerPythia8(Long_t size, void *p);
   static void delete_PileUpMergerPythia8(void *p);
   static void deleteArray_PileUpMergerPythia8(void *p);
   static void destruct_PileUpMergerPythia8(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PileUpMergerPythia8*)
   {
      ::PileUpMergerPythia8 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PileUpMergerPythia8 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PileUpMergerPythia8", ::PileUpMergerPythia8::Class_Version(), "modules/PileUpMergerPythia8.h", 40,
                  typeid(::PileUpMergerPythia8), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PileUpMergerPythia8::Dictionary, isa_proxy, 4,
                  sizeof(::PileUpMergerPythia8) );
      instance.SetNew(&new_PileUpMergerPythia8);
      instance.SetNewArray(&newArray_PileUpMergerPythia8);
      instance.SetDelete(&delete_PileUpMergerPythia8);
      instance.SetDeleteArray(&deleteArray_PileUpMergerPythia8);
      instance.SetDestructor(&destruct_PileUpMergerPythia8);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PileUpMergerPythia8*)
   {
      return GenerateInitInstanceLocal((::PileUpMergerPythia8*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PileUpMergerPythia8::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PileUpMergerPythia8::Class_Name()
{
   return "PileUpMergerPythia8";
}

//______________________________________________________________________________
const char *PileUpMergerPythia8::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PileUpMergerPythia8::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PileUpMergerPythia8::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PileUpMergerPythia8::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PileUpMergerPythia8*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PileUpMergerPythia8::Streamer(TBuffer &R__b)
{
   // Stream an object of class PileUpMergerPythia8.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PileUpMergerPythia8::Class(),this);
   } else {
      R__b.WriteClassBuffer(PileUpMergerPythia8::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PileUpMergerPythia8(void *p) {
      return  p ? new(p) ::PileUpMergerPythia8 : new ::PileUpMergerPythia8;
   }
   static void *newArray_PileUpMergerPythia8(Long_t nElements, void *p) {
      return p ? new(p) ::PileUpMergerPythia8[nElements] : new ::PileUpMergerPythia8[nElements];
   }
   // Wrapper around operator delete
   static void delete_PileUpMergerPythia8(void *p) {
      delete ((::PileUpMergerPythia8*)p);
   }
   static void deleteArray_PileUpMergerPythia8(void *p) {
      delete [] ((::PileUpMergerPythia8*)p);
   }
   static void destruct_PileUpMergerPythia8(void *p) {
      typedef ::PileUpMergerPythia8 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PileUpMergerPythia8

namespace {
  void TriggerDictionaryInitialization_Pythia8Dict_Impl() {
    static const char* headers[] = {
0
    };
    static const char* includePaths[] = {
"external",
"/Users/olmo/programs/root-6.10.08/include",
"/Users/olmo/programs/Delphes-3.4.1/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Pythia8Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$modules/PileUpMergerPythia8.h")))  PileUpMergerPythia8;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Pythia8Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PileUpMergerPythia8.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PileUpMergerPythia8+;

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"PileUpMergerPythia8", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Pythia8Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Pythia8Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Pythia8Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Pythia8Dict() {
  TriggerDictionaryInitialization_Pythia8Dict_Impl();
}
