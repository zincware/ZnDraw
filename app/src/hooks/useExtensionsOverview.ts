import { useQuery } from "@tanstack/react-query";
import {
  getRoomExtensionsOverview,
  getGlobalExtensionsOverview,
  getExtensionDetailedAnalytics,
} from "../myapi/client";

export const useRoomExtensionsOverview = (
  roomId: string,
  filters?: { category?: string; search?: string }
) => {
  return useQuery({
    queryKey: ["extensions", "room", roomId, filters],
    queryFn: () => getRoomExtensionsOverview(roomId, filters),
    staleTime: 30000, // 30 seconds
    refetchInterval: 30000,
  });
};

export const useGlobalExtensionsOverview = (
  filters?: { category?: string; search?: string }
) => {
  return useQuery({
    queryKey: ["extensions", "global", filters],
    queryFn: () => getGlobalExtensionsOverview(filters),
    staleTime: 60000, // 1 minute
  });
};

export const useExtensionDetailedAnalytics = (
  roomId: string,
  category: string,
  extension: string
) => {
  return useQuery({
    queryKey: ["extensions", "analytics", roomId, category, extension],
    queryFn: () => getExtensionDetailedAnalytics(roomId, category, extension),
    staleTime: 60000,
  });
};
