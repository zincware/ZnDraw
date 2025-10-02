Overview

  Change frame index storage from "{physical_index}" to "{source_room_id}:{physical_index}" to enable:
  - Template frames shared across multiple rooms
  - On-demand storage creation
  - Clear frame provenance
  - Reduced memory duplication

  ---
  Phase 1: Validation & Safety

  1.1 Room/Template ID Validation
  - Location: events.py:on_join() (line 199-247)
  - Action: Add validation to reject room IDs containing ":"
  if ":" in room:
      return {"success": False, "error": "Room ID cannot contain ':'"}
  - Also apply to:
    - promote_room_to_template() - validate template IDs
    - Any other room creation entry points

  ---
  Phase 2: Storage Management

  2.1 Update get_storage() Function
  - Location: routes.py:get_storage() (line 42)
  - Current: Already accepts any room_id - no changes needed!
  - Note: Storage is lazily created on first access

  2.2 Storage Lifetime
  - Keep in memory: Templates persist (never cleaned)
  - Consider cleanup: Room storage when room is deleted/empty (future enhancement)

  ---
  Phase 3: Template Initialization

  3.1 Initialize Room with Template Frames
  - Location: events.py:on_join() after line 246
  - Action: When creating a new room with a template:
  # After: r.set(f"room:{room}:template", template_id)

  if template_id != "empty":
      # Copy template's frame mappings to new room
      template_indices_key = f"room:{template_id}:trajectory:indices"
      template_frames = r.zrange(template_indices_key, 0, -1, withscores=True)

      if template_frames:
          room_indices_key = f"room:{room}:trajectory:indices"
          # Add template frames with updated format
          for member, score in template_frames:
              # Parse existing format (might be "old_template:idx" or just "idx")
              if ":" in member:
                  # Already in new format
                  r.zadd(room_indices_key, {member: score})
              else:
                  # Old format - convert to new
                  r.zadd(room_indices_key, {f"{template_id}:{member}": score})

  ---
  Phase 4: Frame Operations

  4.1 Update get_frames() (line 66)
  - Current:
  physical_index = int(frame_mapping[frame_id])
  frame_data = storage.get(physical_index, keys=requested_keys)
  - New:
  mapping_entry = frame_mapping[frame_id].decode() if isinstance(frame_mapping[frame_id], bytes) else frame_mapping[frame_id]

  if ":" in mapping_entry:
      source_room_id, physical_index_str = mapping_entry.split(":", 1)
      physical_index = int(physical_index_str)
      source_storage = get_storage(source_room_id)
  else:
      # Backward compatibility: old format (just index)
      physical_index = int(mapping_entry)
      source_storage = storage  # Use current room's storage

  frame_data = source_storage.get(physical_index, keys=requested_keys)

  4.2 Update append_frame() - All Actions
  - append (line 310-326):
  # Change: r.zadd(indices_key, {str(new_physical_index): logical_position})
  # To:
  r.zadd(indices_key, {f"{room_id}:{new_physical_index}": logical_position})
  - replace (line 217-242):
  # Change: pipeline.zadd(indices_key, {str(new_physical_index): target_frame_id})
  # To:
  pipeline.zadd(indices_key, {f"{room_id}:{new_physical_index}": target_frame_id})
  - extend (line 244-274):
  # Change: new_mapping = {str(start_physical_pos + i): start_logical_pos + i ...}
  # To:
  new_mapping = {f"{room_id}:{start_physical_pos + i}": start_logical_pos + i ...}
  - insert (line 276-308):
  # Change: pipeline.zadd(indices_key, {str(new_physical_index): insert_position})
  # To:
  pipeline.zadd(indices_key, {f"{room_id}:{new_physical_index}": insert_position})

  4.3 Update handle_delete_frame() (line 572)
  - Add validation before deletion:
  # After: physical_index_to_remove = int(frame_mapping[frame_id])
  mapping_entry = frame_mapping[frame_id].decode() if isinstance(frame_mapping[frame_id], bytes) else frame_mapping[frame_id]

  # Check if this is a template frame
  template_id = r.get(f"room:{room}:template")
  if ":" in mapping_entry:
      source_room_id = mapping_entry.split(":", 1)[0]
      if source_room_id != room:
          # This is a template frame or from another room
          return {
              "success": False,
              "error": f"Cannot delete template frame. This frame belongs to '{source_room_id}'",
              "error_type": "PermissionError"
          }

  ---
  Phase 5: Edge Cases & Migration

  5.1 Backward Compatibility
  - Support reading old format (just "0", "1") in get_frames()
  - Gradually migrate by rewriting indices during operations (optional)

  5.2 Template Deletion Protection
  - Already handled: Templates are locked permanently (line 488)
  - Future: Add reference counting if template deletion is needed

  5.3 Empty Template Handling
  - "empty" template has no frames, so no special handling needed

  ---
  Testing Strategy

  1. Unit Tests:
    - Room ID validation (reject ":")
    - Template initialization (frame copying)
    - Mixed frame indices (template + room frames)
    - Delete protection (reject template frame deletion)
  2. Integration Tests:
    - Create room from template → verify frames accessible
    - Append to templated room → verify new frames use room_id:idx
    - Delete template frame → verify rejection
    - Delete room frame → verify success
  3. Migration Test:
    - Existing room with old format → verify still readable
    - Perform operation → verify migrates to new format

  ---
  Potential Issues

  1. Colon in room IDs - Strictly validate to prevent issues
  2. Storage memory - Templates persist in memory (acceptable tradeoff)
  3. Migration complexity - Support both formats during transition
  4. Template deletion - Currently locked permanently (good)
  5. Redis memory - Frame mappings now slightly larger (minimal impact)

  ---
  Summary

  The approach is well-designed and achieves your goals. Key points:

  ✅ Storage efficiency - Templates shared, not duplicated✅ Clear separation - template_id:idx vs room_id:idx✅ Lazy loading - Storage created
  on-demand✅ Backward compatible - Can support old format during migration✅ Safe deletion - Template frames protected

  Recommended order:
  1. Phase 1 (Validation) - Prevent issues early
  2. Phase 3 (Template Init) - Set up frame copying
  3. Phase 4 (Operations) - Update frame handling
  4. Phase 5 (Migration) - Handle legacy data

# Comment
You don't need to be backwards compatible, you MUST ignore all parts that are purely for backwards compatibility and migration from previous data 
formats. Your plan otherwise looks solid, please implement. For the rejection of `:` you can simply replace the `:` with a `_` and update the token in 
the `Client` and log a user warning. 