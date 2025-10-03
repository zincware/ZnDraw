import { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import Box from '@mui/material/Box';
import Container from '@mui/material/Container';
import Typography from '@mui/material/Typography';
import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableContainer from '@mui/material/TableContainer';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';
import Paper from '@mui/material/Paper';
import Button from '@mui/material/Button';
import CircularProgress from '@mui/material/CircularProgress';
import Alert from '@mui/material/Alert';

interface Template {
  id: string;
  name: string;
  description: string;
}

export default function TemplateSelectionPage() {
  const [templates, setTemplates] = useState<Template[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const navigate = useNavigate();

  useEffect(() => {
    const fetchTemplates = async () => {
      try {
        const response = await fetch('/api/templates');
        if (!response.ok) {
          throw new Error(`Failed to fetch templates: ${response.status}`);
        }
        const data = await response.json();

        // Filter out "empty" template for decision logic
        const nonEmptyTemplates = data.filter((t: Template) => t.id !== 'empty');

        if (nonEmptyTemplates.length === 0) {
          // Only "empty" template exists, auto-forward to empty room
          const roomUuid = crypto.randomUUID();
          const userUuid = crypto.randomUUID();
          navigate(`/rooms/${roomUuid}/${userUuid}?template=empty`);
        } else if (nonEmptyTemplates.length === 1) {
          // Only one non-empty template exists, auto-forward to it
          const roomUuid = crypto.randomUUID();
          const userUuid = crypto.randomUUID();
          navigate(`/rooms/${roomUuid}/${userUuid}?template=${nonEmptyTemplates[0].id}`);
        } else {
          // Multiple non-empty templates, show selection table
          setTemplates(data);
        }
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error');
      } finally {
        setLoading(false);
      }
    };

    fetchTemplates();
  }, [navigate]);

  const handleJoinTemplate = (templateId: string) => {
    const roomUuid = crypto.randomUUID();
    const userUuid = crypto.randomUUID();
    navigate(`/rooms/${roomUuid}/${userUuid}?template=${templateId}`);
  };

  if (loading) {
    return (
      <Container maxWidth="md">
        <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '100vh' }}>
          <CircularProgress />
        </Box>
      </Container>
    );
  }

  if (error) {
    return (
      <Container maxWidth="md">
        <Box sx={{ mt: 4 }}>
          <Alert severity="error">{error}</Alert>
        </Box>
      </Container>
    );
  }

  return (
    <Container maxWidth="md">
      <Box sx={{ mt: 4, mb: 4 }}>
        <Typography variant="h3" component="h1" gutterBottom>
          Welcome to Collaborative Viewer
        </Typography>
        <Typography variant="subtitle1" color="text.secondary" gutterBottom>
          Select a template to create a new room
        </Typography>

        <TableContainer component={Paper} sx={{ mt: 3 }}>
          <Table>
            <TableHead>
              <TableRow>
                <TableCell><strong>Name</strong></TableCell>
                <TableCell><strong>Description</strong></TableCell>
                <TableCell align="right"><strong>Action</strong></TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {templates.map((template) => (
                <TableRow key={template.id} hover>
                  <TableCell>{template.name}</TableCell>
                  <TableCell>{template.description}</TableCell>
                  <TableCell align="right">
                    <Button
                      variant="contained"
                      color="primary"
                      onClick={() => handleJoinTemplate(template.id)}
                    >
                      Join
                    </Button>
                  </TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </TableContainer>
      </Box>
    </Container>
  );
}
