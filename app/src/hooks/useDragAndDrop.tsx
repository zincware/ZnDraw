import { useState, useCallback, DragEvent } from "react";
import { uploadFile, getSupportedTypes } from "../myapi/client";
import { useAppStore } from "../store";
import { useQueryClient } from "@tanstack/react-query";

interface UseDragAndDropReturn {
  isDragging: boolean;
  handleDragOver: (e: DragEvent<HTMLDivElement>) => void;
  handleDragEnter: (e: DragEvent<HTMLDivElement>) => void;
  handleDragLeave: (e: DragEvent<HTMLDivElement>) => void;
  handleDrop: (e: DragEvent<HTMLDivElement>) => void;
}

export function useDragAndDrop(): UseDragAndDropReturn {
  const [isDragging, setIsDragging] = useState(false);
  const [dragCounter, setDragCounter] = useState(0);
  // Use individual selectors to prevent unnecessary re-renders
  const setRoomId = useAppStore((state) => state.setRoomId);
  const queryClient = useQueryClient();

  const handleDragOver = useCallback((e: DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
  }, []);

  const handleDragEnter = useCallback((e: DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setDragCounter((prev) => prev + 1);
    setIsDragging(true);
  }, []);

  const handleDragLeave = useCallback((e: DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setDragCounter((prev) => {
      const newCount = prev - 1;
      if (newCount === 0) {
        setIsDragging(false);
      }
      return newCount;
    });
  }, []);

  const isSupportedFormat = useCallback(async (ext: string): Promise<boolean> => {
    try {
      const types = await getSupportedTypes();
      return types.extensions.includes(`.${ext.toLowerCase()}`);
    } catch (error) {
      console.error("Error checking supported formats:", error);
      return false;
    }
  }, []);

  const handleDrop = useCallback(
    async (e: DragEvent<HTMLDivElement>) => {
      e.preventDefault();
      e.stopPropagation();
      setIsDragging(false);
      setDragCounter(0);

      if (!e.dataTransfer.files || e.dataTransfer.files.length === 0) {
        return;
      }

      const files = Array.from(e.dataTransfer.files);

      for (const file of files) {
        const ext = file.name.split(".").pop()?.toLowerCase();
        if (!ext) {
          console.error(`Cannot determine file extension for: ${file.name}`);
          continue;
        }

        const supported = await isSupportedFormat(ext);
        if (!supported) {
          alert(`Unsupported format: .${ext}\n\nPlease use a supported file format (XYZ, PDB, CIF, H5MD, etc.)`);
          continue;
        }

        try {
          const response = await uploadFile({ file });

          // If a room was created, navigate to it
          if (response.room) {
            setRoomId(response.room);
            // Invalidate queries to refresh room list
            queryClient.invalidateQueries({ queryKey: ["rooms"] });

            // Navigate to the new room
            window.location.href = `/rooms/${response.room}`;
          }
        } catch (error: any) {
          console.error("Upload error:", error);
          const errorMsg = error.response?.data?.error || error.message || "Unknown error";
          alert(`Failed to upload ${file.name}: ${errorMsg}`);
        }
      }
    },
    [setRoomId, queryClient, isSupportedFormat]
  );

  return {
    isDragging,
    handleDragOver,
    handleDragEnter,
    handleDragLeave,
    handleDrop,
  };
}
