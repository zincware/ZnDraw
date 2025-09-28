import { createBrowserRouter, RouterProvider } from 'react-router-dom';
import MainPage from './pages/landingPage';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'

const queryClient = new QueryClient()


const LandingPage = () => {
  return (
    <div className="app-layout">
      <header>
        <h1>Welcome to Collaborative Viewer</h1>
        <p>Please select or create a room to join.</p>
        <a href="/room/testroom/testuser">Join Test Room</a>
      </header>
    </div>
  );
};

const router = createBrowserRouter([
  {
    path: '/',
    element: <LandingPage />, // A landing page with instructions or a room list
  },
  {
    path: '/room/:roomId/:userId',
    element: <MainPage />,
  }
]);

export function App() {
  return <QueryClientProvider client={queryClient}><RouterProvider router={router} /></QueryClientProvider>;
}