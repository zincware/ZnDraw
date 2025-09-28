import { useQuery } from '@tanstack/react-query';

const fetchSchemas = async (room: string, option: string) => {
  const response = await fetch(`/api/rooms/${room}/schema/${option}`);
  if (!response.ok) {
    throw new Error(`Failed to fetch schema for ${option}`);
  }
  return response.json();
};

export const useSchemas = (room: string, option: string) => {
  return useQuery({
    queryKey: ['schemas', room, option],
    queryFn: () => fetchSchemas(room, option),
    staleTime: Infinity, // These schemas are static, they never go stale
    enabled: !!room && !!option, // Only run the query if room and option are available
  });
};