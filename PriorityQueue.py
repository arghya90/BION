import numpy as nmp

class PriorityQueue():
    '''
    The arguments passed to a PriorityQueue must consist of
    objects than can be compared using <.
    Use a tuple (priority, item) if necessary.
    --- Taken form HW07 of 8650 Data Structures ---
    '''

    def __init__(self):
        self._array = []

    def push(self, arg_obj):
        # append at end and bubble up
        self._array.append( arg_obj )
        n = len(self._array)
        self._bubble_up(n-1)
        
    def pop(self):
        n = len(self._array)
        if n==0:
            return None
        if n==1:
            return self._array.pop()
        
        # replace with last item and sift down:
        obj = self._array[0]
        self._array[0] = self._array.pop()
        self._sift_down(0)
        return obj
    
    def _parent(self, arg_n):
        return (arg_n-1)//2

    def _left_child(self, arg_n):
        return 2*arg_n + 1

    def _right_child(self, arg_n):
        return 2*arg_n + 2

    def _bubble_up(self, arg_index):
        index = arg_index   # simple copy
        while index>0:
            cur_item = self._array[index]
            parent_idx = self._parent(index)
            parent_item = self._array[parent_idx]
            
            if cur_item < parent_item:
                # swap with parent
                self._array[parent_idx] = cur_item
                self._array[index] = parent_item
                index = parent_idx
            else:
                break
    
    def _sift_down(self,arg_index):
        index = arg_index   # simple copy
        n = len(self._array)
        
        while index<n:           
            cur_item = self._array[index]
            lc = self._left_child(index)
            if n <= lc:
                break

            # first set small child to left child:
            small_child_item = self._array[lc]
            small_child_idx = lc
            
            # right exists and is smaller?
            rc = self._right_child(index)
            if rc < n:
                r_item = self._array[rc]
                if r_item < small_child_item:
                    # right child is smaller than left child:
                    small_child_item = r_item
                    small_child_idx = rc
            
            # done: we are smaller than both children:
            if cur_item <= small_child_item:
                break
            
            # swap with smallest child:
            self._array[index] = small_child_item
            self._array[small_child_idx] = cur_item
            
            # continue with smallest child:
            index = small_child_idx
        
    def size(self):
        return len(self._array)
    
    def is_empty(self):
        return len(self._array) == 0
    

    def heapify(self, arg_items):
        """ Take an array of unsorted items and replace the contents
        of this priority queue by them. """
        # cleaning the present PQ
        self._array.clear()
        
        #fill the array
        for it in arg_items:
            self._array.append(it)
        
        #heapifying the unsorted input
        n = len(self._array)
        
        idx = n-1
        parent_idx = self._parent(idx)
        while ( parent_idx >= 0 ):
            self._sift_down(parent_idx)
            parent_idx -= 1
            
        return
