#ifndef MMMALLOCATOR_H
#define MMMALLOCATOR_H

/**
 * An aligned memory STL allocator, 
 *
 * Modified from the Mallocator from Stephan T. Lavavej, by some Don guy.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 */
 
// The following headers are required for all allocators.
#include <stddef.h>  // Required for size_t and ptrdiff_t and NULL
#include <new>       // Required for placement new and std::bad_alloc
#include <stdexcept> // Required for std::length_error

#if defined(_MSC_VER)
#include <malloc.h>
void *MMMalloc( size_t alignment, size_t size ){ return _aligned_malloc(size, alignment); }
void MMFree( void *buf ){ return _aligned_free(buf); }
#else
#include <stdlib.h>
void * MMMalloc( size_t alignment, size_t size ){ return aligned_alloc(alignment, size); }
void MMFree( void *buf ){ return free(buf); }
#endif // not defined(_MSC_VER)

// The following headers contain stuff that _mm_Mallocator uses.

template <typename T, std::size_t Alignment>
class _mm_Mallocator{
public:

	// The following will be the same for virtually all allocators.
	typedef T * pointer;
	typedef const T * const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T value_type;
	typedef std::size_t size_type;
	typedef ptrdiff_t difference_type;

	T * address(T& r) const {
		return &r;
	}

	const T * address(const T& s) const {
		return &s;
	}

	std::size_t max_size() const {
		// The following has been carefully written to be independent of
		// the definition of size_t and to avoid signed/unsigned warnings.
		return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
	}


	// The following must be the same for all allocators.
	template <typename U>
	struct rebind{
		typedef _mm_Mallocator<U, Alignment> other;
	} ;

	bool operator!=(const _mm_Mallocator& other) const {
		return !(*this == other);
	}

	void construct(T * const p, const T& t) const {
		void * const pv = static_cast<void *>(p);

		new (pv) T(t);
	}

	void destroy(T * const p) const {
		p->~T();
	}

	// Returns true if and only if storage allocated from *this
	// can be deallocated from other, and vice versa.
	// Always returns true for stateless allocators.
	bool operator==(const _mm_Mallocator& other) const {
		return true;
	}


	// Default constructor, copy constructor, rebinding constructor, and destructor.
	// Empty for stateless allocators.
	_mm_Mallocator() { }

	_mm_Mallocator(const _mm_Mallocator&) { }

	template <typename U> _mm_Mallocator(const _mm_Mallocator<U, Alignment>&) { }

	~_mm_Mallocator() { }


	// The following will be different for each allocator.
	T * allocate(const std::size_t n) const {
		// The return value of allocate(0) is unspecified.
		// Mallocator returns NULL in order to avoid depending
		// on malloc(0)'s implementation-defined behavior
		// (the implementation can define malloc(0) to return NULL,
		// in which case the bad_alloc check below would fire).
		// All allocators can return NULL in this case.
		if (n == 0) {
			return NULL;
		}

		// All allocators should contain an integer overflow check.
		// The Standardization Committee recommends that std::length_error
		// be thrown in the case of integer overflow.
		if (n > max_size()) {
			throw std::length_error("_mm_Mallocator<T>::allocate() - Integer overflow.");
		}

		// _mm_Mallocator wraps _mm_malloc().
		// fix for aligned_alloc refusing to service unaligned buf sizes ...
		size_t buf_size = n * sizeof(T); // overflow is guarded against above
		size_t padded_buf_size = (((buf_size-1)/Alignment)+1) * Alignment; // some bit ops for 1<<N alignment
		void * const pv = MMMalloc(Alignment, padded_buf_size);
		// More info about the new aligned new(Align) T[N] and, worse, delete[](Align) T on https://www.cppstories.com/2019/08/newnew-align/

		// Allocators should throw std::bad_alloc in the case of memory allocation failure.
		if (pv == NULL){
			throw std::bad_alloc();
		}

		return static_cast<T *>(pv);
	}

	void deallocate(T * const p, const std::size_t n) const{
		 // _mm_Mallocator wraps _mm_free().
		MMFree(p);
	}

	// The following will be the same for all allocators that ignore hints.
	template <typename U>
	T * allocate(const std::size_t n, const U * /* const hint */) const {
		return allocate(n);
	}

	// Allocators are not required to be assignable, so
	// all allocators should have a private unimplemented
	// assignment operator. Note that this will trigger the
	// off-by-default (enabled under /Wall) warning C4626
	// "assignment operator could not be generated because a
	// base class assignment operator is inaccessible" within
	// the STL headers, but that warning is useless.
private:
	_mm_Mallocator& operator=(const _mm_Mallocator&);
};
#endif
