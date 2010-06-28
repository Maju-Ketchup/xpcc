// coding: utf-8
// ----------------------------------------------------------------------------
/* Copyright (c) 2009, Roboterclub Aachen e.V.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Roboterclub Aachen e.V. nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY ROBOTERCLUB AACHEN E.V. ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ROBOTERCLUB AACHEN E.V. BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * $Id$
 */
// ----------------------------------------------------------------------------

#ifndef	XPCC__DOUBLY_LINKED_LIST_HPP
#define	XPCC__DOUBLY_LINKED_LIST_HPP

#include <stdint.h>

namespace xpcc
{
	/**
	 * \brief	doubly-linked list
	 * 
	 * 
	 * \todo	implementation
	 * 
	 * \tparam	T	type of list entries
	 * 
	 * \author	Fabian Greif
	 * \ingroup	container
	 */
	template <typename T>
	class DoublyLinkedList
	{
	public:
		DoublyLinkedList();
		
		~DoublyLinkedList();
		
		/// check if there are any nodes in the list
		inline bool
		isEmpty() const;
		
		/// Insert in front
		void
		prepend(const T& value);

		/// Insert at the end of the list
		void
		append(const T& value);
		
		/// Remove the first entry
		void
		removeFront();
		
		void
		removeBack();
		
		/**
		 * \return the first node in the list
		 */
		inline const T&
		getFront() const;

		/**
		 * \return the last node in the list
		 */
		inline const T&
		getBack() const;
		
		/**
		 * \brief	Access the node at position \a index
		 * 
		 * \warning	The implementation has to iterate through the list
		 * 			until it reaches the desired position. Therefore an
		 * 			access via iterator is preferred.
		 */
		const T&
		at(int index) const;
		
	protected:
		struct Node
		{
			T data;
			Node *previous;
			Node *next;
		};
		
		
	private:
		DoublyLinkedList(const DoublyLinkedList& other);
		
		DoublyLinkedList&
		operator = (const DoublyLinkedList& other);
	};
}

#include "doubly_linked_list_impl.hpp"

#endif	// XPCC__DOUBLY_LINKED_LIST_HPP