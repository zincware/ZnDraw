import { createBrowserRouter, RouterProvider, Navigate, useParams } from 'react-router-dom';
import MainPage from './pages/landingPage';
import TemplateSelectionPage from './pages/templateSelection';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';

const queryClient = new QueryClient();

const theme = createTheme({
  colorSchemes: {
    dark: true,
  },
  defaultColorScheme: 'light',
});

// Redirect component for /rooms/:roomId -> /rooms/:roomId/:uuid
const RoomRedirect = () => {
  const { roomId } = useParams<{ roomId: string }>();
  const userId = crypto.randomUUID();
  return <Navigate to={`/rooms/${roomId}/${userId}`} replace />;
};

const router = createBrowserRouter([
  {
    path: '/',
    element: <TemplateSelectionPage />,
  },
  {
    path: '/rooms/:roomId',
    element: <RoomRedirect />,
  },
  {
    path: '/rooms/:roomId/:userId',
    element: <MainPage />,
  },
  {
    path: '/room/:roomId/:userId',
    element: <MainPage />,
  }
]);

export function App() {
  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <QueryClientProvider client={queryClient}>
        <RouterProvider router={router} />
      </QueryClientProvider>
    </ThemeProvider>
  );
}