/*******************************************************************************

                               LinkedList.h
                               ------------

  ****************************************************************************
  *                                                                          *
  *   This library is free software; you can redistribute it and/or          *
  *   modify it under the terms of the GNU Lesser General Public             *
  *   License as published by the Free Software Foundation; either           *
  *   version 2.1 of the License, or (at your option) any later version.     *  
  *                                                                          *
  *   This library is distributed in the hope that it will be useful,        *
  *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      *
  *   Lesser General Public License for more details.                        * 
  *                                                                          *
  *   You should have received a copy of the GNU Lesser General Public       * 
  *   License along with this library; if not, write to the Free Software    *
  *   Foundation, Inc.,                                                      *
  *   59 Temple Place, Suite 330, Boston, MA  02111-1307  USA                *
  *                                                                          *
  ****************************************************************************

  Copyright (C) 2005 The University of Chicago

  Authors: 
  Andrés Colubri

  Description: 
  Defines and implements the TNode, TLinkedList_basic and TLinkedList classes, 
  which allow to construct double-linked lists.

*******************************************************************************/

#ifndef __LINKED_LIST_H__
#define __LINKED_LIST_H__

#include <cstddef> // Needed to use NULL.

#include <iostream>
using namespace std;

const int MAX_INT = 2147483647;

/******************************** class TNode **********************************
The node element of the linked list. The ID information can be used to give a
tree structure to the list formed with these nodes.
*******************************************************************************/
template <class T>
class TNode
{
public:
    TNode();
    ~TNode();

    T *Data() { return pData; }
    int id(int level) { return ID[level]; } // Caution! No level-checking: max valid level == IDlength-1
    int *idArray() { return ID; }
    int idLength() { return IDlength; }
    int idMaxLevel() { return IDlength - 1; }

    TNode *Next() { return pNext; }
    TNode *Prev() { return pPrev; }

    void SetData(T *pNewData) { if (pData == NULL) pData = pNewData; else *pData = *pNewData; }
    void ChangeData(T *pNewData) { *pData = *pNewData; }

    void SetID(int id0);
    void SetID(int id0, int id1);
    void SetID(int id0, int id1, int id2);
    void SetID(int *NewID, int NewIDlength);

    void SetNext(TNode *pNewNext) { pNext = pNewNext; }
    void SetPrev(TNode *pNewPrev) { pPrev = pNewPrev; }

    // The overloaded function idCheck returns true if the id of this node matches
    // the first elements specified (in case that the id length is valid).
    bool idCheck(int id0) { return (1 <= IDlength) && (ID[0] == id0); }
    bool idCheck(int id0, int id1) { return (2 <= IDlength) && (ID[0] == id0) && (ID[1] == id1); }
    bool idCheck(int id0, int id1, int id2)
    { return (3 <= IDlength) && (ID[0] == id0) && (ID[1] == id1) && (ID[2] == id2);}
    bool idCheck(int *SomeIDs, int MaxLevel);

    // The overloaded function idCheckAll returns true if the id of this node matches
    // all the elements specified (in case that the id length is valid).
    bool idCheckAll(int id0) { return (1 == IDlength) && (ID[0] == id0); }
    bool idCheckAll(int id0, int id1) { return (2 == IDlength) && (ID[0] == id0) && (ID[1] == id1); }
    bool idCheckAll(int id0, int id1, int id2)
    { return (3 == IDlength) && (ID[0] == id0) && (ID[1] == id1) && (ID[2] == id2);}
    bool idCheckAll(int *AllIDs, int MaxLevel);
private:
    T *pData; // Pointer to node data.

    // This pointer is used to construct a node ID of the form id0.id1.id2... where
    // ID[0] == id0, ID[1] == id1, ID[2] == id2, etc., are integer numbers GREATER OR EQUAL TO ZERO.
    // Each index of the id can be used to define the position of the Node in the level of a tree.
    int *ID;
    int IDlength; // Number of levels in the ID.

    TNode *pNext; // Next node.
    TNode *pPrev; // Previous node.

    void DestroyData() { delete pData; pData = NULL; }
    void DestroyID() { delete[] ID; ID = NULL; IDlength = 0; }
};

// Multiple-index hash table structure.
template <class T>
struct THashTableMultRec
{
    int length; // Number of nodes in the current and lower levels.
    int lengthNNull; // Number of nodes in the current level with non-null data.
    int firstNNullID; // First ID in the current or lower level with non-null data.
    int count; // Number of leaf nodes below this node.
    bool leaf; // Indicates if the position is a node or not.
    int idx; // Data index.
    TNode<T> *pNode; // Pointer to node.
};

/************************ class TLinkedList_basic *******************************
This class stores and manages a double-linked list. Using the multimple-index hash
table, this list can be logically organized and accessed as if it where a tree.
To do this, assign to each node an ID chain of the form id0.id1.id2...
Then call CreateHashTable(). After doing any structure change to the list (not
content), the hash table must be re-created.

Supose that we have 6 nodes and we assing the id's 1.1, 1.2, 2.1, 2.4, 1.1.1 and
2.4.0, then the tree structure of the list will be:

Level 0          /-------------------------\
                 |                         |
Level 1    /------------\             /----------\
           |            |             |          |
       Node0(1.1)   Node1(1.1)    Node2(2.1) Node3(2.4) <-- Node's ids don't need to be consecutive.
           |                                     |
Level 2    |                                     |
           |                                     |
       Node4(1.1.1)                          Node5(2.4.0)

The "real" underlying structure of the list is just Node0<->Node1...<->Node5
*******************************************************************************/
template <class T>
class TLinkedList_basic
{
public:
    TLinkedList_basic(); // Default constructor.
    TLinkedList_basic(const TLinkedList_basic<T> &list); // Copy constructor.
    ~TLinkedList_basic(); // Default destructor.

    // Methods to retrieve data. Important note: the indexing is 0-based.
    int Count() { return NodeCount; } // Number of nodes.
    int CurrentIndex() { return idxCurrent; } // Index of the node currently selected.
    T *CurrentData() // Pointer to the data of the current node.
    { if (pCurrent != NULL) return pCurrent->Data(); else return NULL; }

    // Inserts a new node between the nodes idx-1 and idx. If idx is 0, the node
    // is inserted at the beginning of the list (queued). If idx == NodeCount,
    // then the node is inserted at the end (added or pushed). Returns the index
    // of the new node. pCurrent points to the new node.
    // *Very important note when adding new elements*: To the functions Insert/Add
    // a pointer must be passed, and this pointer must not be deleted outside the list, because
    // an error would be generated (trying to delete the same pointer twice).
    // Avoid passing to Insert/Add the address of a static variable:
    // T data; list.Add(&data);
    // because data is deleted when destroying the list and also when the current
    // scope finishes.
    int Insert(int idx, T *pNewData);
    int Insert(int idx, T *pNewData, int id0)
    { int n = Insert(idx, pNewData); SetCurrentID(id0); return n; }
    int Insert(int idx, T *pNewData, int id0, int id1)
    { int n = Insert(idx, pNewData); SetCurrentID(id0, id1); return n; }
    int Insert(int idx, T *pNewData, int id0, int id1, int id2)
    { int n = Insert(idx, pNewData); SetCurrentID(id0, id1, id2); return n; }
    int Insert(int idx, T *pNewData, int *NewID, int NewIDlength)
    { int n = Insert(idx, pNewData); SetCurrentID(NewID, NewIDlength); return n; }

    // Inserts a new node at the end of the list.
    int Add(T *pNewData) { return Insert(NodeCount, pNewData); }
    int Add(T *pNewData, int id0)
    { int n = Add(pNewData); SetCurrentID(id0); return n; }
    int Add(T *pNewData, int id0, int id1)
    { int n = Add(pNewData); SetCurrentID(id0, id1); return n; }
    int Add(T *pNewData, int id0, int id1, int id2)
    { int n = Add(pNewData); SetCurrentID(id0, id1, id2); return n; }
    int Add(T *pNewData, int *NewID, int NewIDlength)
    { int n = Add(pNewData); SetCurrentID(NewID, NewIDlength); return n; }

    void Delete(int idx); // Deletes the node idx. Doesn't change, if possible, the current node.
    void DeleteRange(int idx0, int idx1); // Deletes [idx0, idx1].
    void DeleteLast() { Delete(NodeCount-1); }

    void Clear(); // Deletes the entire list, but not the hash table.
    void ClearAll(); // Deletes the entire list and the hash table.

    // Change the node data.
    void SetCurrentData(T *pNewData); // The current node.
    void SetData(int idx, T *pNewData) { GoToNode(idx); SetCurrentData(pNewData); } // The node idx.

    // Change the ID of the current node.
    void SetCurrentID(int id0) { if (pCurrent != NULL) pCurrent->SetID(id0); }
    void SetCurrentID(int id0, int id1) { if (pCurrent != NULL) pCurrent->SetID(id0, id1); }
    void SetCurrentID(int id0, int id1, int id2) { if (pCurrent != NULL) pCurrent->SetID(id0, id1, id2); }
    void SetCurrentID(int *NewID, int NewIDlength) { if (pCurrent != NULL) pCurrent->SetID(NewID, NewIDlength); }

    // Change the ID of the node idx.
    void SetID(int idx, int id0) { GoToNode(idx); SetCurrentID(id0); }
    void SetID(int idx, int id0, int id1) { GoToNode(idx); SetCurrentID(id0, id1); }
    void SetID(int idx, int id0, int id1, int id2) { GoToNode(idx); SetCurrentID(id0, id1, id2); }
    void SetID(int idx, int *NewID, int NewIDlength) { GoToNode(idx); SetCurrentID(NewID, NewIDlength); }

    // Give the ID of the current node.
    void GetCurrentID(int &id0) { if (pCurrent != NULL) id0 = pCurrent->id(0); }
    void GetCurrentID(int &id0, int &id1)
    { if (pCurrent != NULL) { id0 = pCurrent->id(0); id1 = pCurrent->id(1); } }
    void GetCurrentID(int &id0, int &id1, int &id2)
    { if (pCurrent != NULL) { id0 = pCurrent->id(0); id1 = pCurrent->id(1); id2 = pCurrent->id(2); } }
    void GetCurrentID(int **ID, int &IDlength)
    { if (pCurrent != NULL) { *ID = pCurrent->idArray(); IDlength = pCurrent->idLength(); } }

    // Give the ID of the node idx.
    void IDXToID(int idx, int &id0) { GoToNode(idx); GetCurrentID(id0); }
    void IDXToID(int idx, int &id0, int &id1) { GoToNode(idx); GetCurrentID(id0, id1);  }
    void IDXToID(int idx, int &id0, int &id1, int &id2) { GoToNode(idx); GetCurrentID(id0, id1, id2);  }
    void IDXToID(int idx, int **ID, int &IDlength) { GoToNode(idx); GetCurrentID(ID, IDlength);  }

    // Check if the id of the current node matches the first elements specified.
    bool EqualSomeIDs(int id0) { if (pCurrent != NULL) return pCurrent->idCheck(id0); else return false; }
    bool EqualSomeIDs(int id0, int id1) { if (pCurrent != NULL) return pCurrent->idCheck(id0, id1); else return false; }
    bool EqualSomeIDs(int id0, int id1, int id2)
    { if (pCurrent != NULL) return pCurrent->idCheck(id0, id1, id2); else return false; }
    bool EqualSomeIDs(int *SomeIDs, int MaxLevel)
    { if (pCurrent != NULL) return pCurrent->idCheck(SomeIDs, MaxLevel); else return false; }

    // Check if the current node has the specified ID.
    bool EqualID(int id0) { if (pCurrent != NULL) return pCurrent->idCheckAll(id0); else return false; }
    bool EqualID(int id0, int id1) { if (pCurrent != NULL) return pCurrent->idCheckAll(id0, id1); else return false; }
    bool EqualID(int id0, int id1, int id2)
    { if (pCurrent != NULL) return pCurrent->idCheckAll(id0, id1, id2); else return false; }
    bool EqualID(int *AllIDs, int MaxLevel)
    { if (pCurrent != NULL) return pCurrent->idCheckAll(AllIDs, MaxLevel); else return false; }

    // Returns the maximum id at level, given that the previous id's are PrevIDs[0]...
    // PrevIDs[level-1]. For instance, if level == 2 and PrevIDs == {1, 2}, then
    // it will give max { x : nodes with id == 1.2.x... }. Doesn't need the hash tables.
    int MaxIDatLevel(int *PrevIDs, int PrevLevel);
    int MaxAllIDLength(); // Gives the maximum ID length calculated over all the nodes.

    // Methods to traverse the list. All of them give the pointer to the node
    // data which results current after the movement, NULL if an invalid position
    // was specified.
    T *GoToNext();
    T *GoToPrev();
    T *GoToFirst();
    T *GoToLast();
    T *GoToNode(int idx);

    void CreateHashTable(); // Constructs the hash table
    void DestroyHashTable(); // Frees the memory used by the hash table.
    bool UsingHashTable() { return (HashTableI != NULL) && (HashTableM != NULL); }

    // Hashing functions.
    int HashFuncI(int n)
    { return n; }
    int HashFuncM(int id0)
    { return id0 + 1; }
    int HashFuncM(int id0, int id1)
    { return id0 * MaxIDLength[1] + id1 + 1 + NElements[0]; }
    int HashFuncM(int id0, int id1, int id2)
    { return MaxIDLength[2] * (id0 * MaxIDLength[1] + id1) + id2 + 1 +
             NElements[0] + NElements[1]; }
    int HashFuncM(int *ID, int IDlength);

    // To access the multiple-index hash table directly.
    int HashTableMSize() { return HTMSize; }
    THashTableMultRec<T> HashTableMRec(int idx) { return HashTableM[idx]; }
    int HashTableMRec_length(int idx) { return HashTableM[idx].length; }
    int HashTableMRec_lengthNNull(int idx) { return HashTableM[idx].lengthNNull; }
    int HashTableMRec_firstNNullID(int idx) { return HashTableM[idx].firstNNullID; }
    int HashTableMRec_count(int idx) { return HashTableM[idx].count; }
    int HashTableMRec_leaf(int idx) { return HashTableM[idx].leaf; }
    int HashTableMRec_idx(int idx) { return HashTableM[idx].idx; }
    TNode<T> *HashTableMRec_pNode(int idx) { return HashTableM[idx].pNode; }

    // Returns the data of node idx using the hash table. No range checking is performed!
    T *GoToNodeI(int idx) { return HashTableI[idx]->Data(); }

    // Check if the id of the node idx matches the first elements specified.
    bool EqualSomeIDsI(int idx, int id0) { return HashTableI[idx]->idCheck(id0); }
    bool EqualSomeIDsI(int idx, int id0, int id1) { return HashTableI[idx]->idCheck(id0, id1); }
    bool EqualSomeIDsI(int idx, int id0, int id1, int id2)
    { return HashTableI[idx]->idCheck(id0, id1, id2); }
    bool EqualSomeIDsI(int idx, int *SomeIDs, int MaxLevel)
    { return HashTableI[idx]->idCheck(SomeIDs, MaxLevel); }

    // Check if the node idx has the specified ID.
    bool EqualIDsI(int idx, int id0) { return HashTableI[idx]->idCheckAll(id0); }
    bool EqualIDsI(int idx, int id0, int id1) { return HashTableI[idx]->idCheckAll(id0, id1); }
    bool EqualIDsI(int idx, int id0, int id1, int id2)
    { return HashTableI[idx]->idCheckAll(id0, id1, id2); }
    bool EqualIDsI(int idx, int *AllIDs, int MaxLevel)
    { return HashTableI[idx]->idCheckAll(AllIDs, MaxLevel); }

    // Give the ID of node idx, using the hash table.
    void IDXToIDsI(int idx, int &id0) { id0 = HashTableI[idx]->id(0); }
    void IDXToIDsI(int idx, int &id0, int &id1)
    { id0 = HashTableI[idx]->id(0); id1 = HashTableI[idx]->id(1); }
    void IDXToIDsI(int idx, int &id0, int &id1, int &id2)
    { id0 = HashTableI[idx]->id(0); id1 = HashTableI[idx]->id(1); id2 = HashTableI[idx]->id(2); }
    void IDXToIDsI(int idx, int **ID, int &IDlength)
    { *ID = HashTableI[idx]->idArray(); IDlength = HashTableI[idx]->idLength(); }

    // Note about the (overloaded) GoToNodeM function: no range checking on the id's is performed!
    // nor is checked whether the hash table was initialized or not.
    T *GoToNodeM(int id0) // Returns data of node id0.
    { return HashTableM[HashFuncM(id0)].pNode->Data(); }
    T *GoToNodeM(int id0, int id1) // Returns data of node id0.id1.
    { return HashTableM[HashFuncM(id0, id1)].pNode->Data(); }
    T *GoToNodeM(int id0, int id1, int id2) // Returns data of node id0.id1.id2.
    { return HashTableM[HashFuncM(id0, id1, id2)].pNode->Data(); }
    T *GoToNodeM(int *ID, int IDlength) // Returns data of node ID[0].ID[1]...ID[IDlength-1]
    { return HashTableM[HashFuncM(ID, IDlength)].pNode->Data(); }

    // Return true if the specified node is a leaf in the logical tree tree (has no descendants).
    bool LeafNode() { return HashTableM == NULL; }
    bool LeafNode(int id0) { return HashTableM[HashFuncM(id0)].leaf; }
    bool LeafNode(int id0, int id1) { return HashTableM[HashFuncM(id0, id1)].leaf; }
    bool LeafNode(int id0, int id1, int id2) { return HashTableM[HashFuncM(id0, id1, id2)].leaf; }
    bool LeafNode(int *ID, int IDlength) { return HashTableM[HashFuncM(ID, IDlength)].leaf; }

    // Give the index of the node with the specified ID. The hash tree is required. To do it
    // without the hash tree, the following method can be used:
    // idx = Count() - 1;
    // GoToNode(idx);
    // while (!EqualID(...) && (-1 < idx)) { idx--; GoToNode(idx); }
    // return idx;
    // -1 is given if no match is found.
    int IDsToIDX(int id0) { return HashTableM[HashFuncM(id0)].idx; }
    int IDsToIDX(int id0, int id1) { return HashTableM[HashFuncM(id0, id1)].idx; }
    int IDsToIDX(int id0, int id1, int id2) { return HashTableM[HashFuncM(id0, id1, id2)].idx; }
    int IDsToIDX(int *ID, int IDlength) { return HashTableM[HashFuncM(ID, IDlength)].idx; }

    // Returns the level length, given that the previous id's are PrevIDs[0]...
    // PrevIDs[level-1]. Note that MaxIDatLevel+1 == IDLevelLength. This function can be used
    // only after creating the hash tree. Caution! no range checking is performed.
    int IDLevelLength() // Level 0.
    { if (HashTableM != NULL) return HashTableM[0].length; else return 0; }
    int IDLevelLength(int id0) // Level 1.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length))
          return HashTableM[HashFuncM(id0)].length; else return 0; }
    int IDLevelLength(int id0, int id1) // Level 2.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length))
          return HashTableM[HashFuncM(id0, id1)].length; else return 0; }
    int IDLevelLength(int id0, int id1, int id2) // Level 3.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length) &
          (id2 < HashTableM[HashFuncM(id0, id1)].length))
          return HashTableM[HashFuncM(id0, id1, id2)].length; else return 0; }
    // Doesn't control that the provided ID's are valid (as in the previous cases).
    int IDLevelLength(int *PrevIDs, int PrevLevel)
    { return HashTableM[HashFuncM(PrevIDs, PrevLevel)].length; }

    // Returns the level non-null length, given that the previous id's are PrevIDs[0]...
    // PrevIDs[level-1]. Note that MaxIDatLevel+1 == IDLevelLength. This function can be used
    // only after creating the hash tree. Caution! no range checking is performed.
    int IDLevelLengthNNull() // Level 0.
    { if (HashTableM != NULL) return HashTableM[0].lengthNNull; else return 0; }
    int IDLevelLengthNNull(int id0) // Level 1.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length))
          return HashTableM[HashFuncM(id0)].lengthNNull; else return 0; }
    int IDLevelLengthNNull(int id0, int id1) // Level 2.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length))
          return HashTableM[HashFuncM(id0, id1)].lengthNNull; else return 0; }
    int IDLevelLengthNNull(int id0, int id1, int id2) // Level 3.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length) &
          (id2 < HashTableM[HashFuncM(id0, id1)].length))
          return HashTableM[HashFuncM(id0, id1, id2)].lengthNNull; else return 0; }
    int IDLevelLengthNNull(int *PrevIDs, int PrevLevel)
    { return HashTableM[HashFuncM(PrevIDs, PrevLevel)].lengthNNull; }

    // Returns the level count, given that the previous id's are PrevIDs[0]...
    // PrevIDs[level-1]. This function can be used only after creating the hash
    // tree. Caution! no range checking is performed.
    int IDLevelCount() // Level 0.
    { if (HashTableM != NULL) return HashTableM[0].count; else return 0; }
    int IDLevelCount(int id0) // Level 1.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length))
          return HashTableM[HashFuncM(id0)].count; else return 0; }
    int IDLevelCount(int id0, int id1) // Level 2.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length))
          return HashTableM[HashFuncM(id0, id1)].count; else return 0; }
    int IDLevelCount(int id0, int id1, int id2) // Level 3.
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length) &
          (id2 < HashTableM[HashFuncM(id0, id1)].length))
          return HashTableM[HashFuncM(id0, id1, id2)].count; else return 0; }
    int IDLevelCount(int *PrevIDs, int PrevLevel)
    { return HashTableM[HashFuncM(PrevIDs, PrevLevel)].count; }

    // Give the first valid id value for each level. Hash tree required.
    int FirstIDAtLevel() // First ID at level 0. Check valid ID's.
    { if ((HashTableM != NULL) && (0 < HashTableM[0].length))
          return HashTableM[0].firstNNullID; else return 0; }
    int FirstIDAtLevel(int id0)
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length))
          return HashTableM[HashFuncM(id0)].firstNNullID; else return 0; }
    int FirstIDAtLevel(int id0, int id1)
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length))
          return HashTableM[HashFuncM(id0, id1)].firstNNullID; else return 0; }
    int FirstIDAtLevel(int id0, int id1, int id2)
    { if ((HashTableM != NULL) && (id0 < HashTableM[0].length) &&
          (id1 < HashTableM[HashFuncM(id0)].length) &&
          (id2 < HashTableM[HashFuncM(id0, id1)].length))
          return HashTableM[HashFuncM(id0, id1, id2)].firstNNullID; else return 0; }
    int FirstIDAtLevel(int *PrevIDs, int PrevLevel) // First ID at level. (doesn't check valid ID's).
    { return HashTableM[HashFuncM(PrevIDs, PrevLevel)].firstNNullID; }

    // Change the node data by accessing the nodes through the hash table.
    void ChangeDataI(int idx, T *pNewData) { TNode<T> *p = HashTableI[idx]; p->ChangeData(pNewData); }

    // Change the node data by accessing the nodes through the hash tree.
    void ChangeDataM(int id0, T *pNewData)
    { THashTableMultRec<T> *p = &HashTableM[HashFuncM(id0)]; p->pNode->ChangeData(pNewData); }
    void ChangeDataM(int id0, int id1, T *pNewData)
    { THashTableMultRec<T> *p = &HashTableM[HashFuncM(id0, id1)]; p->pNode->ChangeData(pNewData); }
    void ChangeDataM(int id0, int id1, int id2, T *pNewData)
    { THashTableMultRec<T> *p = &HashTableM[HashFuncM(id0, id1, id2)]; p->pNode->ChangeData(pNewData); }
    void ChangeDataM(int *ID, int IDlength, T *pNewData)
    { THashTableMultRec<T> *p = &HashTableM[HashFuncM(ID, IDlength)]; p->pNode->ChangeData(pNewData); }

    // Overloaded operators.
    T *operator [] (int idx) { return GoToNodeI(idx); }
    T *operator () (int id0) { return GoToNodeM(id0); }
    T *operator () (int id0, int id1) { return GoToNodeM(id0, id1); }
    T *operator () (int id0, int id1, int id2) { return GoToNodeM(id0, id1, id2); }

    // Unary addition operator: the nodes of list are appended to this list.
    // Read this: <------------------------------------------ !!!!!!!!!!!!!!!!!!
    // If T is a class type, the assignation operator for T should be overloaded
    // to take into account pointer copy, etc.
    TLinkedList_basic<T> &operator += (const TLinkedList_basic<T> &list);
    TLinkedList_basic<T> &operator = (const TLinkedList_basic<T> &list); // Assignement operator.
    // Binary addition operator: the nodes of list are appended to this list and
    // the resulting list is returned.
    TLinkedList_basic<T> operator + (const TLinkedList_basic<T> &list);
protected:
    TNode<T> *pFirst; // First node of the list.
    TNode<T> *pLast; // Last node of the list.
    TNode<T> *pCurrent; // Current node.

    int idxCurrent; // Index of the current node (0-based).
    int NodeCount; // Number of nodes.

    TNode<T> **HashTableI; // Identity-index hash table. Its size is equal to the number of nodes.
    THashTableMultRec<T> *HashTableM; // Multiple-index hash table.
    int HTMSize; // Size of the multiple-index hash table.

    int *MaxIDLength; // Maximum ID-lengths array. Element l contains the maximum length for level l.
    int *NElements; // Contains the number of elements at each level.

    void InitData();
    void InitMaxIDLength(int MaxLength); // Initializes MaxIDLength.
    void InitNElements(int MaxLength); // Initializes NElements.

    void CreateHashTableI(); // Constructs the identity-index hash table.
    void DestroyHashTableI(); // Destroys the identity-index hash table.

    void CreateHashTableM(); // Constructs the multiple-index hash table.
    void DestroyHashTableM(); // Destroys the multiple-index hash table.
};

/*************************** class TLinkedList *********************************
This is TLinkedList_basic "flattened" to avoid dealing with pointers.
*******************************************************************************/
template <class T>
class TLinkedList : public TLinkedList_basic<T>
{
public:
    // Derived constructors.
    TLinkedList() : TLinkedList_basic<T>() {}
    TLinkedList(const TLinkedList<T> &list) : TLinkedList_basic<T>(list) {}
    // New constructors.
    TLinkedList(const T &nulldata) : TLinkedList_basic<T>() { NULLdata = nulldata; }

    T CurrentData() { if (TLinkedList_basic<T>::pCurrent != NULL) return *TLinkedList_basic<T>::pCurrent->Data(); else return NULLdata; }

    int Insert(int idx, const T &NewData) { return TLinkedList_basic<T>::Insert(idx, new T(NewData)); }
    int Insert(int idx, const T &NewData, int id0)
    { int n = Insert(idx, NewData); SetCurrentID(id0); return n; }
    int Insert(int idx, const T &NewData, int id0, int id1)
    { int n = Insert(idx, NewData); SetCurrentID(id0, id1); return n; }
    int Insert(int idx, const T &NewData, int id0, int id1, int id2)
    { int n = Insert(idx, NewData); SetCurrentID(id0, id1, id2); return n; }
    int Insert(int idx, const T &NewData, int *NewID, int NewIDlength)
    { int n = Insert(idx, NewData); SetCurrentID(NewID, NewIDlength); return n; }

    int Add(const T &NewData) { return Insert(TLinkedList_basic<T>::NodeCount, NewData); }
    int Add(const T &NewData, int id0)
    { int n = Add(NewData); SetCurrentID(id0); return n; }
    int Add(const T &NewData, int id0, int id1)
    { int n = Add(NewData); SetCurrentID(id0, id1); return n; }
    int Add(const T &NewData, int id0, int id1, int id2)
    { int n = Add(NewData); SetCurrentID(id0, id1, id2); return n; }
    int Add(const T &NewData, int *NewID, int NewIDlength)
    { int n = Add(NewData); SetCurrentID(NewID, NewIDlength); return n; }

    void SetCurrentData(const T &NewData)
    { if (TLinkedList_basic<T>::pCurrent != NULL) TLinkedList_basic<T>::SetCurrentData(new T(NewData)); }
    void SetData(int idx, const T &NewData) { GoToNode(idx); SetCurrentData(NewData); }

    void SetCurrentID(int id0) { TLinkedList_basic<T>::SetCurrentID(id0); }
    void SetCurrentID(int id0, int id1) { TLinkedList_basic<T>::SetCurrentID(id0, id1); }
    void SetCurrentID(int id0, int id1, int id2) { TLinkedList_basic<T>::SetCurrentID(id0, id1, id2); }
    void SetCurrentID(int *NewID, int NewIDlength) { TLinkedList_basic<T>::SetCurrentID(NewID, NewIDlength); }

    T GoToNext() { return *TLinkedList_basic<T>::GoToNext(); }
    T GoToPrev() { return *TLinkedList_basic<T>::GoToPrev(); }
    T GoToFirst() { return *TLinkedList_basic<T>::GoToFirst(); }
    T GoToLast() { return *TLinkedList_basic<T>::GoToLast(); }
    T GoToNode(int idx) { return *TLinkedList_basic<T>::GoToNode(idx); }

    T GoToNodeI(int idx) { return *TLinkedList_basic<T>::HashTableI[idx]->Data(); }

    T GoToNodeM(int id0) { return *TLinkedList_basic<T>::GoToNodeM(id0); }
    T GoToNodeM(int id0, int id1) { return *TLinkedList_basic<T>::GoToNodeM(id0, id1); }
    T GoToNodeM(int id0, int id1, int id2) { return *TLinkedList_basic<T>::GoToNodeM(id0, id1, id2); }
    T GoToNodeM(int *ID, int IDlength) { return *TLinkedList_basic<T>::GoToNodeM(ID, IDlength); }

    void ChangeDataI(int idx, T &NewData) { TLinkedList_basic<T>::ChangeDataI(idx, &NewData); }

    void ChangeDataM(int id0, T &NewData) { TLinkedList_basic<T>::ChangeDataM(id0, &NewData); }
    void ChangeDataM(int id0, int id1, T &NewData) { TLinkedList_basic<T>::ChangeDataM(id0, id1, &NewData); }
    void ChangeDataM(int id0, int id1, int id2, T &NewData)
    { TLinkedList_basic<T>::ChangeDataM(id0, id1, id2, &NewData); }
    void ChangeDataM(int *ID, int IDlength, T &NewData)
    { TLinkedList_basic<T>::ChangeDataM(ID, IDlength, &NewData); }

    T operator [] (int idx) { return GoToNodeI(idx); }
    T operator () (int id0) { return GoToNodeM(id0); }
    T operator () (int id0, int id1) { return GoToNodeM(id0, id1); }
    T operator () (int id0, int id1, int id2) { return GoToNodeM(id0, id1, id2); }
protected:
    T NULLdata; // Value that is returned when there is no valid data.
};

/********************************* Utilities **********************************/

// Gives true if id is stored in some of the positions of validIDs.
bool ValidID(int id, int *validIDs, int NvalidIDs);

// The implementations cannot be in a separate file because the use of templates.

/**************************** TNode implementation ****************************/

template <class T>
TNode<T>::TNode()
{
    ID = NULL;
    IDlength = 0;
    pData = NULL;
    pNext = NULL;
    pPrev = NULL;
}

template <class T>
TNode<T>::~TNode()
{
    DestroyID();
    DestroyData();
    pNext = NULL;
    pPrev = NULL;
}

template <class T>
void TNode<T>::SetID(int id0)
{
    if (IDlength != 1)
    {
        DestroyID();
        IDlength = 1;
        ID = new int[IDlength];
    }
    ID[0] = id0;
}

template <class T>
void TNode<T>::SetID(int id0, int id1)
{
    if (IDlength != 2)
    {
        DestroyID();
        IDlength = 2;
        ID = new int[IDlength];
    }
    ID[0] = id0; ID[1] = id1;
}

template <class T>
void TNode<T>::SetID(int id0, int id1, int id2)
{
    if (IDlength != 3)
    {
        DestroyID();
        IDlength = 3;
        ID = new int[IDlength];
    }
    ID[0] = id0; ID[1] = id1; ID[2] = id2;
}

template <class T>
void TNode<T>::SetID(int *NewID, int NewIDlength)
{
    if (IDlength != NewIDlength)
    {
        DestroyID();
        ID = NewID;
        IDlength = NewIDlength;
    }
    else
    {
        int i;
        for (i = 0; i < IDlength; i++) ID[i] = NewID[i];
    }
}

template <class T>
bool TNode<T>::idCheck(int *SomeIDs, int MaxLevel)
{
    if (MaxLevel < IDlength)
    {
        int i;
        bool b = true;

        for (i = 0; i <= MaxLevel; i++)
            if (ID[i] != SomeIDs[i])
            {
                b = false;
                break;
            }

        return b;
    }
    else return false;
}

template <class T>
bool TNode<T>::idCheckAll(int *AllIDs, int MaxLevel)
{
    if (MaxLevel == IDlength - 1)
    {
        int i;
        bool b = true;

        for (i = 0; i <= MaxLevel; i++)
            if (ID[i] != AllIDs[i])
            {
                b = false;
                break;
            }

        return b;
    }
    else return false;
}

/********************* TLinkedList_basic implementation ************************/

template <class T>
TLinkedList_basic<T>::TLinkedList_basic() { InitData(); }

template <class T>
TLinkedList_basic<T>::TLinkedList_basic(const TLinkedList_basic<T> &list)
{ InitData(); *this = list; }

template <class T>
TLinkedList_basic<T>::~TLinkedList_basic() { ClearAll(); }

template <class T>
int TLinkedList_basic<T>::Insert(int idx, T *pNewData)
{
    TNode<T> *n0, *n1;

    if ((NodeCount < MAX_INT) && (0 <= idx) && (idx <= NodeCount))
    {

        if (idx == 0)
        {
            // A new element is inserted at the beginning of the list.
            pCurrent = new TNode<T>;
            if (pFirst != NULL)
            {
                pFirst->SetPrev(pCurrent);
                pCurrent->SetNext(pFirst);
            }
            pFirst = pCurrent;
        }
        else
        {
            // A new node is inserted between idx-1 and idx.
            GoToNode(idx - 1);

            n0 = pCurrent; // node idx-1
            n1 = pCurrent->Next(); // node idx

            pCurrent = new TNode<T>;
            pCurrent->SetPrev(n0);
            pCurrent->SetNext(n1);

            n0->SetNext(pCurrent);
            if (n1 != NULL) n1->SetPrev(pCurrent);
        }

        NodeCount++;
        pCurrent->SetData(pNewData);
        if (pCurrent->Next() == NULL) pLast = pCurrent;
        idxCurrent = idx;
        return idx;
    }
    else return -1;
}

template <class T>
void TLinkedList_basic<T>::Delete(int idx)
{
    TNode<T> *TempNode, *Node0;
    int idx0;

    if ((0 <= idx) && (idx <= NodeCount - 1)) // Checking validity of index.
    {
        Node0 = pCurrent;
        idx0 = idxCurrent;
        if (idx < idx0) idx0--;

        if (idx == 0)
        {
            // The base node is deleted (this is done in a separate case because
            // the base node doesn't have a previous node).
            TempNode = pFirst;
            pFirst = pFirst->Next();
            if (pFirst != NULL) pFirst->SetPrev(NULL);
        }
        else
        {
            // The node previous to the one to be deleted is located.
            GoToNode(idx - 1);
            TempNode = pCurrent->Next();
            // TempNode now points to the node with position idx.
            // The node with position idx is removed from the list.
            pCurrent->SetNext(pCurrent->Next()->Next());
            if (pCurrent->Next() != NULL) pCurrent->Next()->SetPrev(pCurrent);
        }

        if (TempNode == Node0)
        {
            // The original position was deleted. Going to the first node.
            pCurrent = pFirst;
            if (pFirst == NULL) idxCurrent = -1; else idxCurrent = 0;
        }
        else
        {
            // Restoring original position.
            pCurrent = Node0;
            idxCurrent = idx0;
        }

        if (TempNode == pLast) pLast = TempNode->Prev();
        NodeCount--;
        delete TempNode;
    }
}

template <class T>
void TLinkedList_basic<T>::DeleteRange(int idx0, int idx1)
{
    for (int idx = idx1; idx0 <= idx; idx--) Delete(idx);
}

template <class T>
void TLinkedList_basic<T>::Clear()
{
    TNode<T> *TempNode;

    while (pFirst != NULL)
    {
        TempNode = pFirst;
        pFirst = pFirst->Next();
        delete TempNode;
    }

    pLast = NULL;
    pCurrent = NULL;
    idxCurrent = -1;
    NodeCount = 0;
}

template <class T>
void TLinkedList_basic<T>::ClearAll()
{
    DestroyHashTable();
    Clear();
}

template <class T>
void TLinkedList_basic<T>::SetCurrentData(T *pNewData)
{
    if (pCurrent != NULL) pCurrent->SetData(pNewData);
}

template <class T>
int TLinkedList_basic<T>::MaxIDatLevel(int *PrevIDs, int PrevLevel)
{
    int max = -1; // MaxID == -1 means that the specified level is empty.
    TNode<T> *p = pFirst;

    while (p != NULL)
    {
        if ((PrevLevel < p->idLength()) &&
            ((PrevLevel - 1 < 0) || p->idCheck(PrevIDs, PrevLevel - 1)) &&
            (max < p->id(PrevLevel))) max = p->id(PrevLevel);
        p = p->Next();
    }

    return max;
}

template <class T>
int TLinkedList_basic<T>::MaxAllIDLength()
{
    int max = 0;
    TNode<T> *p = pFirst;

    while (p != NULL)
    {
        if (max < p->idLength()) max = p->idLength();
        p = p->Next();
    }

    return max;
}

template <class T>
T *TLinkedList_basic<T>::GoToNext()
{
    if ((pCurrent != NULL) && (pCurrent->Next() != NULL))
    {
        pCurrent = pCurrent->Next();
        return pCurrent->Data();
        idxCurrent++;
    }
    else return NULL;
}

template <class T>
T *TLinkedList_basic<T>::GoToPrev()
{
    if ((pCurrent != NULL) && (pCurrent->Prev() != NULL))
    {
        pCurrent = pCurrent->Prev();
        return pCurrent->Data();
        idxCurrent--;
    }
    else return NULL;
}

template <class T>
T *TLinkedList_basic<T>::GoToFirst()
{
    if (pFirst != NULL)
    {
        pCurrent = pFirst;
        return pCurrent->Data();
    }
    else return NULL;
}

template <class T>
T *TLinkedList_basic<T>::GoToLast()
{
    if (pLast != NULL)
    {
        pCurrent = pLast;
        return pCurrent->Data();
    }
    else return NULL;
}

template <class T>
T *TLinkedList_basic<T>::GoToNode(int idx)
{
    int i, dif;

    if ((0 <= idx) && (idx <= NodeCount - 1)) // Checking validity of index.
    {
        dif = idx - idxCurrent;
        if (0 <= dif)
            for (i = 1; i <= dif; i++) pCurrent = pCurrent->Next();
        else
            for (i = -1; dif <= i; i--) pCurrent = pCurrent->Prev();
        idxCurrent = idx;
        return pCurrent->Data();
    }
    else return NULL;
}

template <class T>
void TLinkedList_basic<T>::CreateHashTable()
{
    CreateHashTableI();
    CreateHashTableM();
}

template <class T>
void TLinkedList_basic<T>::DestroyHashTable()
{
    DestroyHashTableI();
    DestroyHashTableM();
}

template <class T>
int TLinkedList_basic<T>::HashFuncM(int *ID, int IDlength)
{
    if (0 < IDlength)
    {
        int n = ID[IDlength - 1] + 1;
        int i, p = MaxIDLength[IDlength - 1];
        for (i = IDlength - 2; 0 <= i; i--)
        {
            n = n + ID[i] * p + NElements[i];
            p *= MaxIDLength[i];
        }
        return n;
    }
    else return 0;
}

template <class T>
TLinkedList_basic<T> &TLinkedList_basic<T>::operator += (const TLinkedList_basic<T> &list)
{
    T *pNewData;

    TNode<T> *p = list.pFirst;
    while (p != NULL)
    {
        pNewData = new T; // Create a new pointer to T for each element of list to add.

        // Read this: <------------------------------------------ !!!!!!!!!!!!!!!!!!
        // If T is a class type, the assignation operator should be overloaded to
        // take into account pointer copy, etc.
        *pNewData = *p->Data();

        Add(pNewData);
        p = p->Next();
    }

   if (UsingHashTable()) CreateHashTable();
   return *this;
}

template <class T>
TLinkedList_basic<T> &TLinkedList_basic<T>::operator = (const TLinkedList_basic<T> &list)
{
    Clear();
    *this += list;
    return *this;
}

template <class T>
TLinkedList_basic<T> TLinkedList_basic<T>::operator + (const TLinkedList_basic<T> &list)
{
   TLinkedList_basic<T> newlist(*this);
   newlist += list;
   return newlist;
}

/* -------------------- Implementation of protected functions --------------- */

template <class T>
void TLinkedList_basic<T>::InitData()
{
    pFirst = NULL;
    pLast = NULL;
    pCurrent = NULL;

    idxCurrent = -1;
    NodeCount = 0;

    HashTableI = NULL;
    HashTableM = NULL;
    HTMSize = 0;

    MaxIDLength = NULL;
    NElements = NULL;
}

template <class T>
void TLinkedList_basic<T>::InitMaxIDLength(int MaxLength)
{
    int i;
    for (i = 0; i < MaxLength; i++) MaxIDLength[i] = -1;

    TNode<T> *p = pFirst;
    while (p != NULL)
    {
        for (i = 0; i < p->idLength(); i++)
            if (MaxIDLength[i] < p->id(i) + 1) MaxIDLength[i] = p->id(i) + 1;
        p = p->Next();
    }
}

template <class T>
void TLinkedList_basic<T>::InitNElements(int MaxLength)
{
    int i;
    NElements[0] = MaxIDLength[0];
    for (i = 1; i < MaxLength; i++) NElements[i] = NElements[i - 1] * MaxIDLength[i];
}

template <class T>
void TLinkedList_basic<T>::CreateHashTableI()
{
    DestroyHashTableI();
    HashTableI = new TNode<T>*[NodeCount];
    int i = 0;
    TNode<T> *p = pFirst;
    while (p != NULL)
    {
        HashTableI[i] = p;
        p = p->Next();
        i++;
    }
}

template <class T>
void TLinkedList_basic<T>::DestroyHashTableI()
{
    if (HashTableI != NULL)
    {
        delete[] HashTableI;
        HashTableI = NULL;
    }
}

template <class T>
void TLinkedList_basic<T>::CreateHashTableM()
{
    DestroyHashTableM();
    int MaxLength = MaxAllIDLength();
    if (0 < MaxLength)
    {
        int i;
        MaxIDLength = new int[MaxLength];
        NElements = new int[MaxLength];
        InitMaxIDLength(MaxLength);
        InitNElements(MaxLength);

        // Size of the hash table.
        HTMSize = 1; // + 1 to store the root node.
        for (i = 0; i < MaxLength; i++) HTMSize += NElements[i];
        HashTableM = new THashTableMultRec<T>[HTMSize];

        // Initializing the hash table.
        for (i = 0; i < HTMSize; i++)
        {
            HashTableM[i].length = 0;
            HashTableM[i].lengthNNull = 0;
            HashTableM[i].firstNNullID = -1;
            HashTableM[i].count = 0;
            HashTableM[i].pNode = NULL;
        }

        // Storing data in the hash table.
        int l, pid, pid1, n = -1;
        TNode<T> *p = pFirst;
        while (p != NULL) // Iterating through nodes.
        {
            n++;
            i = HashFuncM(p->idArray(), p->idLength());
            HashTableM[i].idx = n;
            HashTableM[i].pNode = p;

            i = HashFuncM(p->idArray(), p->idLength() - 1);

            // Updating non-null length.
            HashTableM[i].lengthNNull++;

            // Updating lengths, counts and first non-null ID's and non-leaf elements.
            for (l = 0; l < p->idLength(); l++)
            {
                pid = p->id(l); pid1 = pid + 1;
                i = HashFuncM(p->idArray(), l);

                if (HashTableM[i].length < pid1) HashTableM[i].length = pid1;

                HashTableM[i].count++;

                if ((HashTableM[i].firstNNullID == -1) ||
                    (pid < HashTableM[i].firstNNullID)) HashTableM[i].firstNNullID = pid;
            }

            p = p->Next();
        }

        // Setting firstNNullID, leaf and idx.
        for (i = 0; i < HTMSize; i++)
        {
            if (HashTableM[i].firstNNullID == -1)
            {
                HashTableM[i].firstNNullID = 0;
                HashTableM[i].leaf = true;
            }
            else HashTableM[i].leaf = false;
            if (HashTableM[i].pNode == NULL) HashTableM[i].idx = -1;
        }
    }
}

template <class T>
void TLinkedList_basic<T>::DestroyHashTableM()
{
    if (MaxIDLength != NULL) { delete[] MaxIDLength; MaxIDLength = NULL; }
    if (NElements != NULL) { delete[] NElements; NElements = NULL; }
    if (HashTableM != NULL) { delete[] HashTableM; HashTableM = NULL; }
}

#endif
