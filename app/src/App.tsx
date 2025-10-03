import { createBrowserRouter, RouterProvider } from 'react-router-dom';
import MainPage from './pages/landingPage';
import TemplateSelectionPage from './pages/templateSelection';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'

const queryClient = new QueryClient()

const router = createBrowserRouter([
  {
    path: '/',
    element: <TemplateSelectionPage />,
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
  return <QueryClientProvider client={queryClient}><RouterProvider router={router} /></QueryClientProvider>;
}