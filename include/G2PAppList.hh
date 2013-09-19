// -*- C++ -*-

/* class G2PAppList
 * A Collection of G2PAppBase classes.
 */

// History:
//   Sep 2013, C. Gu, First public version
//

#ifndef G2P_APPLIST_H
#define	G2P_APPLIST_H

#include "TList.h"

class G2PAppBase;

class G2PAppList : public TList {
public:
    G2PAppList();
    virtual ~G2PAppList();

    G2PAppBase* Find(const char* name) const;
    G2PAppList* FindList(const char* name) const;
    G2PAppList* FindList(int priority) const;

private:
    ClassDef(G2PAppList, 2)
};

#endif

