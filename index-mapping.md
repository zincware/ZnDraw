# Improve the Frame Mapping

Each file fetching operations is an $O(N)$ and needs refactoring.

Use A More Scalable Indexing Strategy (The "Gap" Method)!

- stop using contiguous integer scores (0, 1, 2, ...). Instead, use floating-point scores.

**Append**: To append a frame, get the score of the last element and add 1.0 (or a large integer like 1000).

**Insert**: To insert a frame at logical position i, get the scores of the frames at i-1 and i. The new score will be the average of those two scores. Example: to insert between score 2.0 and 3.0, the new score is 2.5. This is a single $O(\log N)$ ZADD operation.

**Delete**: A delete is just a ZREM. No re-numbering is needed. The ordering is preserved.

Pitfalls:

- Singleton lists with only one or no elements need special handling
```py
if not prev_element and not next_element:
    new_score = 1.0  # first frame
elif not prev_element:
    new_score = next_element[0][1] - 1.0
elif not next_element:
    new_score = prev_element[0][1] + 1.0
else:
    new_score = (prev_element[0][1] + next_element[0][1]) / 2.0
```

- Handling Precision Drift
Write a renormalization function that can be called via a /api/rooms/<room_id>/renormalize endpoint. This function will reassign scores to be contiguous integers starting from 0. This is an $O(N)$ operation but should be infrequent.

# Implementation

Consider extracting repeatable logic into reusable functions or classes to enhance maintainability and reduce code duplication.
Create new helper or utility modules if necessary to encapsulate related functionalities.
