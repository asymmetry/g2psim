// -*- C++ -*-

/* class G2PAppList
 * A Collection of G2PAppBase classes.
 */

// History:
//   Sep 2013, C. Gu, First public version
//

#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TClass.h"
#include "TList.h"

#include "G2PAppBase.hh"

#include "G2PAppList.hh"

G2PAppList::G2PAppList()
{
    // Nothing to do
}

G2PAppList::~G2PAppList()
{
    Clear();
}

G2PAppBase *G2PAppList::Find(const char *name) const
{
    static const char *const g2papp = "G2PAppBase";

    TIter next(this);

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        if (aobj->IsA()->InheritsFrom(g2papp) && aobj->IsA()->InheritsFrom(name)) {
            if (!aobj->IsZombie())
                return aobj;
        }
    }

    return NULL;
}

G2PAppList *G2PAppList::FindList(const char *name) const
{
    static const char *const g2papp = "G2PAppBase";

    G2PAppList *list = new G2PAppList();
    TIter next(this);
    int n = 0;

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        if (aobj->IsA()->InheritsFrom(g2papp) && aobj->IsA()->InheritsFrom(name)) {
            if (!aobj->IsZombie()) {
                list->Add(aobj);
                n++;
            }
        }
    }

    return list;
}

G2PAppList *G2PAppList::FindList(int priority) const
{
    static const char *const g2papp = "G2PAppBase";

    G2PAppList *list = new G2PAppList();
    TIter next(this);
    int n = 0;

    while (G2PAppBase *aobj = static_cast<G2PAppBase *>(next())) {
        if (aobj->IsA()->InheritsFrom(g2papp)) {
            if (!aobj->IsZombie()) {
                if (aobj->GetPriority() == priority) {
                    list->Add(aobj);
                    n++;
                }
            }
        }
    }

    return list;
}

ClassImp(G2PAppList)
