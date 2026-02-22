export class MinHeap {
  constructor(maxSize, compareFn = (a, b) => a.key - b.key) {
    this._maxSize = maxSize;
    this._compareFn = compareFn;
    this._heap = [];
  }

  get size() {
    return this._heap.length;
  }

  push(item) {
    if (this._heap.length < this._maxSize) {
      this._heap.push(item);
      this._siftUp(this._heap.length - 1);
    } else if (this._heap.length > 0 && this._compareFn(item, this._heap[0]) < 0) {
      this._heap[0] = item;
      this._siftDown(0);
    }
  }

  toArray() {
    return this._heap.slice().sort(this._compareFn);
  }

  // Internal max-heap: parent is always >= children (by compareFn)
  _siftUp(i) {
    while (i > 0) {
      const parent = (i - 1) >> 1;
      if (this._compareFn(this._heap[i], this._heap[parent]) > 0) {
        [this._heap[i], this._heap[parent]] = [this._heap[parent], this._heap[i]];
        i = parent;
      } else {
        break;
      }
    }
  }

  _siftDown(i) {
    const n = this._heap.length;
    while (true) {
      let largest = i;
      const left = 2 * i + 1;
      const right = 2 * i + 2;
      if (left < n && this._compareFn(this._heap[left], this._heap[largest]) > 0) {
        largest = left;
      }
      if (right < n && this._compareFn(this._heap[right], this._heap[largest]) > 0) {
        largest = right;
      }
      if (largest !== i) {
        [this._heap[i], this._heap[largest]] = [this._heap[largest], this._heap[i]];
        i = largest;
      } else {
        break;
      }
    }
  }
}
